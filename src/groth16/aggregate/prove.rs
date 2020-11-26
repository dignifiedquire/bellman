use std::marker::PhantomData;

use digest::Digest;
use ff::Field;
use groupy::{CurveAffine, CurveProjective};
use itertools::Itertools;

use super::{inner_product, poly::DensePolynomial, structured_scalar_power, SRS};
use crate::bls::Engine;
use crate::groth16::Proof;

pub struct AggregateProof<E: Engine, D: Digest> {
    pub com_a: E::Fqk,
    pub com_b: E::Fqk,
    pub com_c: E::Fqk,
    pub ip_ab: E::Fqk,
    pub agg_c: E::G1,
    pub tipa_proof_ab: PairingInnerProductABProof<E, D>,
    pub tipa_proof_c: MultiExpInnerProductCProof<E, D>,
}

pub struct PairingInnerProductABProof<E: Engine, D: Digest> {
    pub gipa_proof: GIPAProof<E, D>,
    pub final_ck: (E::G2, E::G1), // Key
    pub final_ck_proof: (E::G2, E::G1),
    pub _marker: PhantomData<D>,
}

pub struct GIPAProof<E: Engine, D: Digest> {
    // Contains all the v and w values computed during the recursion
    // where initial w = [g^alpha^2i] and v = [h^beta^2i]
    pub r_commitment_steps: Vec<((E::Fqk, E::Fqk, Vec<E::Fqk>), (E::Fqk, E::Fqk, Vec<E::Fqk>))>, // Output
    // Final A and B value out of the recursion
    pub r_base: (E::G1, E::G2), // Message
    pub _marker: PhantomData<D>,
}
pub struct GIPAAux<E: Engine, D: Digest> {
    // Contains all challenges derived in the recursion
    pub r_transcript: Vec<E::Fr>,
    // Contains final value of v and w computed during recursion
    pub ck_base: (E::G2, E::G1),
    pub _marker: PhantomData<D>,
}

pub struct MultiExpInnerProductCProof<E: Engine, D: Digest> {
    pub gipa_proof: GIPAProofWithSSM<E, D>,
    pub final_ck: E::G2,
    pub final_ck_proof: E::G2,
    pub _marker: PhantomData<D>,
}

pub struct GIPAProofWithSSM<E: Engine, D: Digest> {
    pub r_commitment_steps: Vec<((E::Fqk, Vec<E::G1>), (E::Fqk, Vec<E::G1>))>, // Output
    pub r_base: (E::G1, E::Fr),                                                // Message
    pub _marker: PhantomData<D>,
}

pub struct GIPAAuxWithSSM<E: Engine, D: Digest> {
    pub r_transcript: Vec<E::Fr>,
    pub ck_base: E::G2,
    pub _marker: PhantomData<D>,
}

pub fn aggregate_proofs<E: Engine, D: Digest>(
    ip_srs: &SRS<E>,
    proofs: &[Proof<E>],
) -> AggregateProof<E, D> {
    println!("aggregate proofs");
    let a = proofs
        .iter()
        .map(|proof| proof.a.into_projective())
        .collect::<Vec<E::G1>>();
    let b = proofs
        .iter()
        .map(|proof| proof.b.into_projective())
        .collect::<Vec<E::G2>>();
    let c = proofs
        .iter()
        .map(|proof| proof.c.into_projective())
        .collect::<Vec<E::G1>>();

    // ck_1 = (g^alpha^2*0, ... , g^alpha^2*(n-1)) = v
    // ck_2 = (h^beta^2*0, ... , h^beta^2*(n-1)) = w
    let (ck_1, ck_2) = ip_srs.get_commitment_keys(); // (v, w)

    let com_a = inner_product::pairing::<E>(&a, &ck_1); // T
    let com_b = inner_product::pairing::<E>(&ck_2, &b); // U
    let com_c = inner_product::pairing::<E>(&c, &ck_1); // V

    // Random linear combination of proofs
    let mut counter_nonce: usize = 0;
    let r = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        hash_input.extend_from_slice(&com_a.as_bytes());
        hash_input.extend_from_slice(&com_b.as_bytes());
        hash_input.extend_from_slice(&com_c.as_bytes());
        if let Some(r) =
            E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES])
        {
            break r;
        };

        counter_nonce += 1;
    };

    // Just to easily eyeball values of `r`.
    // let r = {
    //     let mut x = E::Fr::one();
    //     x.add_assign(&E::Fr::one());
    //     x
    // }; // two

    // r_vec = (1, r^2, r^4, ..., r^2n-2)
    let r_vec = structured_scalar_power(proofs.len(), &r);
    println!("r {:?} ", &r);

    // a_r = (A_0^1, A_1^r^2, ..., A_n^r^(2n-2)) = A'
    let a_r = a
        .iter()
        .zip(&r_vec)
        .map(|(a, r)| {
            let mut a = *a;
            a.mul_assign(*r);
            a
        })
        .collect::<Vec<E::G1>>();
    // Z = MUL e( A_i^r_i, B^i )
    let ip_ab = inner_product::pairing::<E>(&a_r, &b);
    // C = MUL C_i ^ r^2i
    let agg_c = inner_product::multiexponentiation::<E::G1>(&c, r_vec.as_slice());

    // ck_1_r = [ v_i^(1/r) ] = [ g^alpha^(2i / r) ]
    let ck_1_r = ck_1
        .iter()
        .zip(&r_vec)
        .map(|(ck, r)| {
            let mut ck = *ck;
            ck.mul_assign(r.inverse().unwrap());
            ck
        })
        .collect::<Vec<E::G2>>(); // v^(r^-1)

    // T = A * v = A^r * v^(r^-1)
    assert_eq!(com_a, inner_product::pairing::<E>(&a_r, &ck_1_r));

    println!("prove with srs shift");
    // TIPP proves that Z = A' * B <=>  a_r * b = ip_ab
    // XXX We dont give Z here maybe optim to avoid recomputing ?
    let tipa_proof_ab = prove_with_srs_shift::<E, D>(&ip_srs, (&a_r, &b), (&ck_1_r, &ck_2), &r);

    println!("prove with structured scalar messages");
    // MIPP
    let tipa_proof_c = prove_with_structured_scalar_message::<E, D>(&ip_srs, (&c, &r_vec), &ck_1);

    AggregateProof {
        com_a, // T
        com_b, // U
        com_c, // V
        ip_ab, // Z
        agg_c,
        tipa_proof_ab,
        tipa_proof_c,
    }
}

/// values is a tuple (A' = A^r, B) and ck is a tuple (v' = v^r-1, w) -
/// it returns a proof that A' * B = Z
// Shifts KZG proof for left message by scalar r (used for efficient composition with aggregation protocols)
// LMC commitment key should already be shifted before being passed as input
fn prove_with_srs_shift<E: Engine, D: Digest>(
    srs: &SRS<E>,
    values: (&[E::G1], &[E::G2]),
    ck: (&[E::G2], &[E::G1]),
    r_shift: &E::Fr,
) -> PairingInnerProductABProof<E, D> {
    // Run GIPA to prove Z = A' * B with adapted commitment keys v' and w
    let (proof, aux) = GIPAProof::<E, D>::prove_with_aux(values, (ck.0, ck.1));

    // Prove final commitment keys are wellformed
    let (ck_a_final, ck_b_final) = aux.ck_base; // final v and w
    let transcript = aux.r_transcript; // all challenges, vector x in the paper
    let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();
    let r_inverse = r_shift.inverse().unwrap(); // r^-1

    // KZG challenge point
    let mut counter_nonce: usize = 0;
    let c = loop {
        // challence  z in paper
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        hash_input.extend_from_slice(&transcript.first().unwrap().as_bytes());
        hash_input.extend_from_slice(ck_a_final.into_affine().into_uncompressed().as_ref());
        hash_input.extend_from_slice(ck_b_final.into_affine().into_uncompressed().as_ref());
        if let Some(c) =
            E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES])
        {
            break c;
        };
        counter_nonce += 1;
    };

    // Complete KZG proofs
    let ck_a_kzg_opening = // opening for v value
        prove_commitment_key_kzg_opening(&srs.h_beta_powers, &transcript_inverse, &r_inverse, &c);
    let ck_b_kzg_opening = // opening for w value - r = 1 in that case since it's not in the formula
        prove_commitment_key_kzg_opening(&srs.g_alpha_powers, &transcript, &<E::Fr>::one(), &c);

    PairingInnerProductABProof {
        gipa_proof: proof,
        final_ck: (ck_a_final, ck_b_final),
        final_ck_proof: (ck_a_kzg_opening, ck_b_kzg_opening),
        _marker: PhantomData,
    }
}

// IP: PairingInnerProduct<E>
// LMC: AFGHOCommitmentG1<E>
// RMC: AFGHOCommitmentG2<E>
// IPC: IdentityCommitment<E::Fqk, E::Fr>
impl<E: Engine, D: Digest> GIPAProof<E, D> {
    pub fn prove_with_aux(
        values: (&[E::G1], &[E::G2]),
        ck: (&[E::G2], &[E::G1]),
    ) -> (Self, GIPAAux<E, D>) {
        let (m_a, m_b) = values;
        let (ck_a, ck_b) = ck;
        Self::_prove((m_a.to_vec(), m_b.to_vec()), (ck_a.to_vec(), ck_b.to_vec()))
    }

    /// Returns vector of recursive commitments and transcripts in reverse order.
    fn _prove(
        values: (Vec<E::G1>, Vec<E::G2>),
        ck: (Vec<E::G2>, Vec<E::G1>),
    ) -> (Self, GIPAAux<E, D>) {
        let (mut m_a, mut m_b) = values;
        let (mut ck_a, mut ck_b) = ck;
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        assert!(m_a.len().is_power_of_two());
        let (m_base, ck_base) = 'recurse: loop {
            if m_a.len() == 1 {
                // base case
                break 'recurse (
                    (m_a[0].clone(), m_b[0].clone()),
                    (ck_a[0].clone(), ck_b[0].clone()),
                );
            } else {
                // recursive step
                // Recurse with problem of half size
                let split = m_a.len() / 2;

                let m_a_right = &m_a[split..];
                let m_a_left = &m_a[..split];
                let ck_a_left = &ck_a[..split];
                let ck_a_right = &ck_a[split..];

                let m_b_left = &m_b[..split];
                let m_b_right = &m_b[split..];
                let ck_b_right = &ck_b[split..];
                let ck_b_left = &ck_b[..split];

                // compute left Tl Ul Zl and right commitments Tr Ur Zr
                let com_left = (
                    // T_left = A [m':] * v [m':] = A_right * v_right
                    inner_product::pairing::<E>(m_a_right, ck_a_left),
                    //  U_left = w [m':] * B [:m'] = w_right * b_left
                    inner_product::pairing::<E>(ck_b_right, m_b_left),
                    // Z_left = A [m':] * B [:m'] = A_right * B_left
                    vec![inner_product::pairing::<E>(m_a_right, m_b_left)],
                );
                let com_right = (
                    // T_right = A [:m'] * v [:m'] = A_left * v_left
                    inner_product::pairing::<E>(m_a_left, ck_a_right),
                    // U_right = w [:m'] * B [m':] = w_left * b_right
                    inner_product::pairing::<E>(ck_b_left, m_b_right),
                    // Z_right = A [:m'] * B [m'] = A_left * B_right
                    vec![inner_product::pairing::<E>(m_a_left, m_b_right)],
                );

                // Fiat-Shamir challenge
                let mut counter_nonce: usize = 0;
                let default_transcript = E::Fr::zero();
                let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (c, c_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    hash_input.extend_from_slice(&transcript.as_bytes());
                    hash_input.extend_from_slice(&com_left.0.as_bytes());
                    hash_input.extend_from_slice(&com_left.1.as_bytes());
                    for c in &com_left.2 {
                        hash_input.extend_from_slice(&c.as_bytes());
                    }
                    hash_input.extend_from_slice(&com_right.0.as_bytes());
                    hash_input.extend_from_slice(&com_right.1.as_bytes());
                    for c in &com_right.2 {
                        hash_input.extend_from_slice(&c.as_bytes());
                    }
                    let c = E::Fr::from_bytes(
                        &D::digest(&hash_input).as_slice()[0..E::Fr::SERIALIZED_BYTES],
                    );
                    if let Some(c) = c {
                        if let Some(c_inv) = c.inverse() {
                            // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                            // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                            break 'challenge (c_inv, c);
                        }
                    }

                    counter_nonce += 1;
                };

                // Set up values for next step of recursion
                // A' = A_right ^ c  ° A_left
                m_a = m_a_right
                    .iter()
                    .map(|a| mul!(*a, c))
                    .zip(m_a_left)
                    .map(|(a_1, a_2)| add!(a_1, a_2))
                    .collect::<Vec<_>>();

                // B' = B_right ^ c^-1  °   B_left
                m_b = m_b_right
                    .iter()
                    .map(|b| mul!(*b, c_inv))
                    .zip(m_b_left)
                    .map(|(b_1, b_2)| add!(b_1, b_2))
                    .collect::<Vec<_>>();

                // v' = v_right^ x^-1  ° v_left
                ck_a = ck_a_right
                    .iter()
                    .map(|a| mul!(*a, c_inv))
                    .zip(ck_a_left)
                    .map(|(a_1, a_2)| add!(a_1, a_2))
                    .collect::<Vec<_>>();

                // w' = w_right^x * w_left
                ck_b = ck_b_right
                    .iter()
                    .map(|b| mul!(*b, c))
                    .zip(ck_b_left)
                    .map(|(b_1, b_2)| add!(b_1, b_2))
                    .collect::<Vec<_>>();

                r_commitment_steps.push((com_left, com_right));
                r_transcript.push(c);
            }
        };
        r_transcript.reverse();
        r_commitment_steps.reverse();
        (
            GIPAProof {
                r_commitment_steps,
                r_base: m_base,
                _marker: PhantomData,
            },
            GIPAAux {
                r_transcript,
                ck_base,
                _marker: PhantomData,
            },
        )
    }
}

// IPAWithSSMProof<
//    MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
//    AFGHOCommitmentG1<P>,
//    IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
//    P,
//    D,
// >;
//
// GIPAProof<
//   IP = MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
//   LMC = AFGHOCommitmentG1<P>,
//   RMC = SSMPlaceholderCommitment<LMC::Scalar>
//   IPC = IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,,
// D>,
//
// IP: MultiexponentiationInnerProduct<<P as PairingEngine>::G1Projective>,
// LMC: AFGHOCommitmentG1<P>,
// RMC: SSMPlaceholderCommitment<LMC::Scalar>,
// IPC: IdentityCommitment<<P as PairingEngine>::G1Projective, <P as PairingEngine>::Fr>,
impl<E: Engine, D: Digest> GIPAProofWithSSM<E, D> {
    pub fn prove_with_aux(
        values: (&[E::G1], &[E::Fr]),
        ck: (&[E::G2], &[E::Fr]),
    ) -> (Self, GIPAAuxWithSSM<E, D>) {
        let (m_a, m_b) = values;
        let (ck_a, ck_t) = ck;
        Self::_prove((m_a.to_vec(), m_b.to_vec()), (ck_a.to_vec(), ck_t.to_vec()))
    }

    /// Returns vector of recursive commitments and transcripts in reverse order.
    fn _prove(
        values: (Vec<E::G1>, Vec<E::Fr>),
        ck: (Vec<E::G2>, Vec<E::Fr>),
    ) -> (Self, GIPAAuxWithSSM<E, D>) {
        let (mut m_a, mut m_b) = values;
        let (mut ck_a, _ck_t) = ck;
        let mut r_commitment_steps = Vec::new();
        let mut r_transcript = Vec::new();
        assert!(m_a.len().is_power_of_two());
        let (m_base, ck_base) = 'recurse: loop {
            if m_a.len() == 1 {
                // base case
                break 'recurse ((m_a[0].clone(), m_b[0].clone()), ck_a[0].clone());
            } else {
                // recursive step
                // Recurse with problem of half size
                let split = m_a.len() / 2;

                let m_a_1 = &m_a[split..];
                let m_a_2 = &m_a[..split];
                let ck_a_1 = &ck_a[..split];
                let ck_a_2 = &ck_a[split..];

                let m_b_1 = &m_b[..split];
                let m_b_2 = &m_b[split..];

                let com_1 = (
                    inner_product::pairing::<E>(m_a_1, ck_a_1),
                    vec![inner_product::multiexponentiation::<E::G1>(m_a_1, m_b_1)],
                );
                let com_2 = (
                    inner_product::pairing::<E>(m_a_2, ck_a_2),
                    vec![inner_product::multiexponentiation::<E::G1>(m_a_2, m_b_2)],
                );

                // Fiat-Shamir challenge
                let mut counter_nonce: usize = 0;
                let default_transcript = E::Fr::zero();
                let transcript = r_transcript.last().unwrap_or(&default_transcript);
                let (c, c_inv) = 'challenge: loop {
                    let mut hash_input = Vec::new();
                    hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
                    hash_input.extend_from_slice(&transcript.as_bytes());
                    hash_input.extend_from_slice(&com_1.0.as_bytes());
                    for c in &com_1.1 {
                        hash_input.extend_from_slice(c.into_affine().into_uncompressed().as_ref());
                    }
                    hash_input.extend_from_slice(&com_2.0.as_bytes());
                    for c in &com_2.1 {
                        hash_input.extend_from_slice(c.into_affine().into_uncompressed().as_ref());
                    }
                    let c = E::Fr::from_bytes(
                        &D::digest(&hash_input).as_slice()[0..E::Fr::SERIALIZED_BYTES],
                    );
                    if let Some(c) = c {
                        if let Some(c_inv) = c.inverse() {
                            // Optimization for multiexponentiation to rescale G2 elements with 128-bit challenge
                            // Swap 'c' and 'c_inv' since can't control bit size of c_inv
                            break 'challenge (c_inv, c);
                        }
                    }

                    counter_nonce += 1;
                };

                // Set up values for next step of recursion
                m_a = m_a_1
                    .iter()
                    .map(|a| mul!(*a, c))
                    .zip(m_a_2)
                    .map(|(a_1, a_2)| add!(a_1, a_2))
                    .collect::<Vec<_>>();

                m_b = m_b_2
                    .iter()
                    .map(|b| mul!(*b, &c_inv))
                    .zip(m_b_1)
                    .map(|(b_1, b_2)| add!(b_1, b_2))
                    .collect::<Vec<_>>();

                ck_a = ck_a_2
                    .iter()
                    .map(|a| mul!(*a, c_inv))
                    .zip(ck_a_1)
                    .map(|(a_1, a_2)| add!(a_1, a_2))
                    .collect::<Vec<_>>();

                r_commitment_steps.push((com_1, com_2));
                r_transcript.push(c);
            }
        };
        r_transcript.reverse();
        r_commitment_steps.reverse();
        (
            GIPAProofWithSSM {
                r_commitment_steps,
                r_base: m_base,
                _marker: PhantomData,
            },
            GIPAAuxWithSSM {
                r_transcript,
                ck_base,
                _marker: PhantomData,
            },
        )
    }
}
/// prove_commitment_key_kzg_opening implements the protocol described in
/// section 6.2.1 which prooves the validity of the final v and w commtiments
/// keys computed during a TIPA proof
pub fn prove_commitment_key_kzg_opening<G: CurveProjective>(
    srs_powers: &Vec<G>,         // [p^e^2i] - can be g^alpha^2i or h^beta^2i
    transcript: &Vec<G::Scalar>, // challenges x from last to first, ^-1 or not
    r_shift: &G::Scalar,         // public randomness used in tipa, ^-1 or not
    kzg_challenge: &G::Scalar,   // specific challenge to this proof
) -> G {
    // f(x) = MUL(j:0 -> l) ( x_(l-j) + r * X^2^(j+1) where l is length of transcript
    let ck_polynomial =
        DensePolynomial::from_coeffs(polynomial_coefficients_from_transcript(transcript, r_shift));
    assert_eq!(srs_powers.len(), ck_polynomial.coeffs().len());
    // TODO: check that g^ck_polynomial.eval(alpha) = w

    // evaluation f(z)
    let ck_polynomial_c_eval =
        polynomial_evaluation_product_form_from_transcript(&transcript, kzg_challenge, &r_shift);

    let mut neg_kzg_challenge = *kzg_challenge;
    neg_kzg_challenge.negate();

    // Q(x)  = f(x) - f(z) / (X - z)
    let quotient_polynomial = &(&ck_polynomial
        - &DensePolynomial::from_coeffs(vec![ck_polynomial_c_eval]))
        / &(DensePolynomial::from_coeffs(vec![neg_kzg_challenge, G::Scalar::one()]));

    let mut quotient_polynomial_coeffs = quotient_polynomial.into_coeffs();
    quotient_polynomial_coeffs.resize(srs_powers.len(), G::Scalar::zero());

    // evaluation of Q(x) at alpha or beta depending on the srs_powers given
    let opening = inner_product::multiexponentiation(srs_powers, &quotient_polynomial_coeffs);
    opening
}

pub(super) fn polynomial_evaluation_product_form_from_transcript<F: Field>(
    transcript: &Vec<F>,
    z: &F,
    r_shift: &F,
) -> F {
    let mut power_2_zr = *z;
    power_2_zr.mul_assign(z);
    power_2_zr.mul_assign(r_shift);
    let mut product_form = Vec::new();
    for x in transcript.iter() {
        product_form.push(add!(F::one(), &mul!(*x, &power_2_zr)));
        power_2_zr.mul_assign(&power_2_zr.clone());
    }
    product_form[1..]
        .iter()
        .fold(product_form[0], |mut acc, curr| {
            acc.mul_assign(curr);
            acc
        })
}

/// Returns coefficients of polynomial in 6.2.1
/// transcript is set of challenges generated: x vector in the paper
/// r_shift is the randomness used for TIPA proof, inversed or not
/// XXX How is that algo is equivalent to the formula in paper ?
fn polynomial_coefficients_from_transcript<F: Field>(transcript: &Vec<F>, r_shift: &F) -> Vec<F> {
    let mut coefficients = vec![F::one()];
    let mut power_2_r = r_shift.clone();
    for (i, x) in transcript.iter().enumerate() {
        // from j in [0..2^i[
        for j in 0..(2_usize).pow(i as u32) {
            // coeffs[j] * (x * r^i)
            coefficients.push(mul!(coefficients[j], &mul!(*x, &power_2_r)));
        }
        // r^2, r^4, r^8...
        power_2_r.mul_assign(&power_2_r.clone());
    }
    // Interleave with 0 coefficients
    coefficients
        .iter()
        .interleave(vec![F::zero()].iter().cycle().take(coefficients.len() - 1))
        .cloned()
        .collect()
}

fn prove_with_structured_scalar_message<E: Engine, D: Digest>(
    srs: &SRS<E>,
    values: (&[E::G1], &[E::Fr]),
    ck: &[E::G2],
) -> MultiExpInnerProductCProof<E, D> {
    // Run GIPA to prove
    let (proof, aux) = GIPAProofWithSSM::<E, D>::prove_with_aux(values, (ck, &vec![])); // TODO: add plaeholder value

    // Prove final commitment key is wellformed
    let ck_a_final = aux.ck_base;
    let transcript = aux.r_transcript;
    let transcript_inverse = transcript.iter().map(|x| x.inverse().unwrap()).collect();

    // KZG challenge point
    let mut counter_nonce: usize = 0;
    let c = loop {
        let mut hash_input = Vec::new();
        hash_input.extend_from_slice(&counter_nonce.to_be_bytes()[..]);
        hash_input.extend_from_slice(&transcript.first().unwrap().as_bytes());
        hash_input.extend_from_slice(ck_a_final.into_affine().into_uncompressed().as_ref());
        if let Some(c) =
            E::Fr::from_bytes(&D::digest(&hash_input).as_slice()[..E::Fr::SERIALIZED_BYTES])
        {
            break c;
        };
        counter_nonce += 1;
    };

    // Complete KZG proof
    let ck_a_kzg_opening = prove_commitment_key_kzg_opening(
        &srs.h_beta_powers,
        &transcript_inverse,
        &E::Fr::one(),
        &c,
    );

    MultiExpInnerProductCProof {
        gipa_proof: proof,
        final_ck: ck_a_final,
        final_ck_proof: ck_a_kzg_opening,
        _marker: PhantomData,
    }
}
