use ff::Field;

macro_rules! mul {
    ($a:expr, $b:expr) => {{
        let mut a = $a;
        a.mul_assign($b);
        a
    }};
}

macro_rules! add {
    ($a:expr, $b:expr) => {{
        let mut a = $a;
        a.add_assign($b);
        a
    }};
}

macro_rules! sub {
    ($a:expr, $b:expr) => {{
        let mut a = $a;
        a.sub_assign($b);
        a
    }};
}

mod inner_product;
mod msm;
mod poly;
mod prove;
mod srs;
mod verify;

pub use self::prove::*;
pub use self::srs::*;
pub use self::verify::*;

/// Returns vector of size num: [1, s^2, s^4, ..., s^(2num-2) ]
fn structured_scalar_power<F: Field>(num: usize, s: &F) -> Vec<F> {
    println!("structured scalar power");
    let mut powers = vec![F::one()];
    let mut s2 = s.clone();
    s2.mul_assign(&s);

    for i in 1..num {
        let mut x = powers[i - 1];
        x.mul_assign(&s2);
        println!("Position {} -> push {:?}", i, &x);
        powers.push(x);
    }
    dbg!(&powers);
    powers
}
