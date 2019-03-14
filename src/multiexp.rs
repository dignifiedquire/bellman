use bit_vec::{self, BitVec};
use pairing::{CurveAffine, CurveProjective, Engine, Field, PrimeField, PrimeFieldRepr};
use rayon;
use std::io;
use std::iter;
use std::sync::Arc;

use super::SynthesisError;

/// An object that builds a source of bases.
pub trait SourceBuilder<G: CurveAffine>: Send + Sync + 'static + Clone {
    type Source: Source<G>;

    fn new(self) -> Self::Source;
}

/// A source of bases, like an iterator.
pub trait Source<G: CurveAffine> {
    /// Parses the element from the source. Fails if the point is at infinity.
    fn add_assign_mixed(
        &mut self,
        to: &mut <G as CurveAffine>::Projective,
    ) -> Result<(), SynthesisError>;

    /// Skips `amt` elements from the source, avoiding deserialization.
    fn skip(&mut self, amt: usize) -> Result<(), SynthesisError>;
}

impl<G: CurveAffine> SourceBuilder<G> for (Arc<Vec<G>>, usize) {
    type Source = (Arc<Vec<G>>, usize);

    fn new(self) -> (Arc<Vec<G>>, usize) {
        (self.0.clone(), self.1)
    }
}

impl<G: CurveAffine> Source<G> for (Arc<Vec<G>>, usize) {
    fn add_assign_mixed(
        &mut self,
        to: &mut <G as CurveAffine>::Projective,
    ) -> Result<(), SynthesisError> {
        if self.0.len() <= self.1 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "expected more bases from source",
            )
            .into());
        }

        if self.0[self.1].is_zero() {
            return Err(SynthesisError::UnexpectedIdentity);
        }

        to.add_assign_mixed(&self.0[self.1]);

        self.1 += 1;

        Ok(())
    }

    fn skip(&mut self, amt: usize) -> Result<(), SynthesisError> {
        if self.0.len() <= self.1 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "expected more bases from source",
            )
            .into());
        }

        self.1 += amt;

        Ok(())
    }
}

pub trait QueryDensity {
    /// Returns whether the base exists.
    type Iter: Iterator<Item = bool>;

    fn iter(self) -> Self::Iter;
    fn get_query_size(self) -> Option<usize>;
}

#[derive(Clone)]
pub struct FullDensity;

impl AsRef<FullDensity> for FullDensity {
    fn as_ref(&self) -> &FullDensity {
        self
    }
}

impl<'a> QueryDensity for &'a FullDensity {
    type Iter = iter::Repeat<bool>;

    fn iter(self) -> Self::Iter {
        iter::repeat(true)
    }

    fn get_query_size(self) -> Option<usize> {
        None
    }
}

pub struct DensityTracker {
    bv: BitVec,
    total_density: usize,
}

impl<'a> QueryDensity for &'a DensityTracker {
    type Iter = bit_vec::Iter<'a>;

    fn iter(self) -> Self::Iter {
        self.bv.iter()
    }

    fn get_query_size(self) -> Option<usize> {
        Some(self.bv.len())
    }
}

impl DensityTracker {
    pub fn new() -> DensityTracker {
        DensityTracker {
            bv: BitVec::new(),
            total_density: 0,
        }
    }

    pub fn add_element(&mut self) {
        self.bv.push(false);
    }

    pub fn inc(&mut self, idx: usize) {
        if !self.bv.get(idx).unwrap() {
            self.bv.set(idx, true);
            self.total_density += 1;
        }
    }

    pub fn get_total_density(&self) -> usize {
        self.total_density
    }
}

fn multiexp_inner<Q, D, G, S>(
    bases: S,
    density_map: D,
    exponents: &[<<G::Engine as Engine>::Fr as PrimeField>::Repr],
    mut skip: u32,
    c: u32,
    handle_trivial: bool,
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G>,
{
    let oldskip = skip;
    skip += c;

    let bases = bases.clone();
    let density_map = density_map.clone();

    let (first, second): (
        Result<<G as CurveAffine>::Projective, SynthesisError>,
        Result<Option<_>, _>,
    ) = rayon::join(
        || {
            let skip = oldskip;
            // Perform this region of the multiexp
            let bases = bases.clone();
            let density_map = density_map.clone();

            // Accumulate the result
            let mut acc = G::Projective::zero();

            // Build a source for the bases
            let mut bases = bases.new();

            // Create space for the buckets
            let mut buckets = vec![<G as CurveAffine>::Projective::zero(); (1 << c) - 1];

            let zero = <G::Engine as Engine>::Fr::zero().into_repr();
            let one = <G::Engine as Engine>::Fr::one().into_repr();

            // Sort the bases into buckets
            for (&exp, density) in exponents.iter().zip(density_map.as_ref().iter()) {
                if density {
                    if exp == zero {
                        bases.skip(1)?;
                    } else if exp == one {
                        if handle_trivial {
                            bases.add_assign_mixed(&mut acc)?;
                        } else {
                            bases.skip(1)?;
                        }
                    } else {
                        let mut exp = exp;
                        exp.shr(skip);
                        let exp = exp.as_ref()[0] % (1 << c);

                        if exp != 0 {
                            bases.add_assign_mixed(&mut buckets[(exp - 1) as usize])?;
                        } else {
                            bases.skip(1)?;
                        }
                    }
                }
            }

            // Summation by parts
            // e.g. 3a + 2b + 1c = a +
            //                    (a) + b +
            //                    ((a) + b) + c
            let mut running_sum = G::Projective::zero();
            for exp in buckets.iter().rev() {
                running_sum.add_assign(&exp);
                acc.add_assign(&running_sum);
            }

            Ok(acc)
        },
        || {
            if skip >= <G::Engine as Engine>::Fr::NUM_BITS {
                // There isn't another region.
                Ok(None)
            } else {
                // There's another region more significant. Calculate and join it with
                // this region recursively.
                multiexp_inner(
                    bases.clone(),
                    density_map.clone(),
                    exponents,
                    skip,
                    c,
                    false,
                )
                .map(|mut higher| {
                    for _ in 0..c {
                        higher.double();
                    }
                    Some(higher)
                })
            }
        },
    );

    match second? {
        Some(mut second) => {
            let first = first?;
            second.add_assign(&first);
            Ok(second)
        }
        None => first,
    }
}

/// Perform multi-exponentiation. The caller is responsible for ensuring the
/// query size is the same as the number of exponents.
pub fn multiexp<Q, D, G, S>(
    bases: S,
    density_map: D,
    exponents: &[<<G::Engine as Engine>::Fr as PrimeField>::Repr],
) -> Result<<G as CurveAffine>::Projective, SynthesisError>
where
    for<'a> &'a Q: QueryDensity,
    D: Send + Sync + 'static + Clone + AsRef<Q>,
    G: CurveAffine,
    S: SourceBuilder<G>,
{
    let c = if exponents.len() < 32 {
        3u32
    } else {
        (f64::from(exponents.len() as u32)).ln().ceil() as u32
    };

    if let Some(query_size) = density_map.as_ref().get_query_size() {
        // If the density map has a known query size, it should not be
        // inconsistent with the number of exponents.

        assert!(query_size == exponents.len());
    }

    multiexp_inner(bases, density_map, exponents, 0, c, true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use pairing::bls12_381::Bls12;
    use pairing::{CurveAffine, CurveProjective, Engine, PrimeField};
    use rand::{self, Rand};
    use rayon::prelude::*;
    use std::sync::Arc;
    use test::{black_box, Bencher};

    const SAMPLES: usize = 1 << 8;

    #[cfg(feature = "profile")]
    use gperftools::profiler::PROFILER;

    #[cfg(feature = "profile")]
    #[inline(always)]
    fn start_profile(stage: &str) {
        PROFILER
            .lock()
            .unwrap()
            .start(format!("./{}.profile", stage))
            .unwrap();
    }

    #[cfg(not(feature = "profile"))]
    #[inline(always)]
    fn start_profile(_stage: &str) {}

    #[cfg(feature = "profile")]
    #[inline(always)]
    fn stop_profile() {
        PROFILER.lock().unwrap().stop().unwrap();
    }

    #[cfg(not(feature = "profile"))]
    #[inline(always)]
    fn stop_profile() {}

    fn naive_multiexp<G: CurveAffine>(
        bases: Arc<Vec<G>>,
        exponents: Arc<Vec<<G::Scalar as PrimeField>::Repr>>,
    ) -> G::Projective {
        assert_eq!(bases.len(), exponents.len());

        let mut acc = G::Projective::zero();

        for (base, exp) in bases.iter().zip(exponents.iter()) {
            acc.add_assign(&base.mul(*exp));
        }

        acc
    }

    fn multiexp_par_iter<G: CurveAffine>(
        bases: Vec<G>,
        exponents: Vec<<G::Scalar as PrimeField>::Repr>,
    ) -> G::Projective {
        assert_eq!(bases.len(), exponents.len());

        bases
            .into_par_iter()
            .zip(exponents.par_iter())
            .map(|(base, exp)| base.mul(*exp))
            .reduce(G::Projective::zero, |mut left, right| {
                left.add_assign(&right);
                left
            })
    }

    #[test]
    fn test_with_bls12() {
        let rng = &mut rand::thread_rng();
        const SAMPLES: usize = 1 << 6; //14;
        let v = Arc::new(
            (0..SAMPLES)
                .map(|_| <Bls12 as Engine>::Fr::rand(rng).into_repr())
                .collect::<Vec<_>>(),
        );
        let g = Arc::new(
            (0..SAMPLES)
                .map(|_| <Bls12 as Engine>::G1::rand(rng).into_affine())
                .collect::<Vec<_>>(),
        );

        let naive = naive_multiexp(g.clone(), v.clone());

        let fast = multiexp((g.clone(), 0), FullDensity, &v).unwrap();

        let par_iter = multiexp_par_iter(g.to_vec(), v.to_vec());

        assert_eq!(naive, fast);
        assert_eq!(naive, par_iter);
    }

    #[bench]
    fn bench_multiexp_naive(b: &mut Bencher) {
        let rng = &mut rand::thread_rng();
        let v = Arc::new(
            (0..SAMPLES)
                .map(|_| <Bls12 as Engine>::Fr::rand(rng).into_repr())
                .collect::<Vec<_>>(),
        );
        let g = Arc::new(
            (0..SAMPLES)
                .map(|_| <Bls12 as Engine>::G1::rand(rng).into_affine())
                .collect::<Vec<_>>(),
        );

        b.iter(|| black_box(naive_multiexp(g.clone(), v.clone())));
    }

    #[bench]
    fn bench_multiexp_pool(b: &mut Bencher) {
        let rng = &mut rand::thread_rng();
        let v = (0..SAMPLES)
            .map(|_| <Bls12 as Engine>::Fr::rand(rng).into_repr())
            .collect::<Vec<_>>();

        let g = Arc::new(
            (0..SAMPLES)
                .map(|_| <Bls12 as Engine>::G1::rand(rng).into_affine())
                .collect::<Vec<_>>(),
        );

        start_profile("pool-multiexp");
        b.iter(|| black_box(multiexp((g.clone(), 0), FullDensity, &v.clone()).unwrap()));
        stop_profile();
    }

    #[bench]
    fn bench_multiexp_par_iter(b: &mut Bencher) {
        let rng = &mut rand::thread_rng();
        let v = (0..SAMPLES)
            .map(|_| <Bls12 as Engine>::Fr::rand(rng).into_repr())
            .collect::<Vec<_>>();
        let g = (0..SAMPLES)
            .map(|_| <Bls12 as Engine>::G1::rand(rng).into_affine())
            .collect::<Vec<_>>();

        start_profile("par-iter");
        b.iter(|| black_box(multiexp_par_iter(g.clone(), v.clone())));
        stop_profile();
    }
}
