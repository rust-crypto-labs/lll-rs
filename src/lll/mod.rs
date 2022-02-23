//! The Lenstra-Lenstra-Lovasz algorithm [LLL82]

use crate::algebra::{BigNum, Float, FromExt, Matrix, Scalar};

/// Lattice reduction using the original Lenstra-Lenstra-Lovasz algorithm
///
/// This implementation uses generic Scalars for arithmetic operations.
/// The value of `delta` is set to 0.75.
///
///   - `basis`: A generating matrix for the lattice
///
/// The basis is reduced in-place.
fn lattice_reduce<S: Scalar>(basis: &mut Matrix<S::Integer>) {
    // Parameter delta in the Lovasz condition
    let delta = S::Fraction::from_ext((3, 4));

    let (n, _) = basis.dimensions();
    let mut swap_condition = true;

    while swap_condition {
        // Perform rounded Gram-Schmidt orthogonalisation
        for i in 0..n {
            for k in 1..i {
                let j = i - k;

                let b_i = &basis[i];
                let b_j = &basis[j];
                let alpha = S::round_div(b_i.dot(b_j), b_j.dot(b_j));
                basis[i] = b_i.sub(&b_j.mulf(alpha));
            }
        }

        // Check for the Lovasz condition and swap columns if appropriate
        swap_condition = false;
        for i in 0..n - 1 {
            let b_i = &basis[i];
            let b_ip1 = &basis[i + 1];

            let lhs = S::Fraction::from_ext(&b_i.dot(b_i)) * &delta;

            let alpha = S::round_div(b_ip1.dot(b_i), b_i.dot(b_i));
            let vec_rhs = b_ip1.add(&b_i.mulf(alpha));
            let rhs = vec_rhs.dot(&vec_rhs);

            if lhs > rhs {
                basis.swap(i, i + 1);
                swap_condition = true;
                break;
            }
        }
    }
}

/// Lattice reduction using the original Lenstra-Lenstra-Lovasz algorithm
///
/// This implementation uses generic `rug::Integer` and `rug::Fraction` for arithmetic operations.
/// The value of `delta` is set to 0.75.
///
///   - `basis`: A generating matrix for the lattice
///
/// The basis is reduced in-place.
#[deprecated(
    note = "Current implementation might yield incorrect results. Use l2.lll_bignum() instead"
)]
pub fn lll_bignum(basis: &mut Matrix<rug::Integer>) {
    lattice_reduce::<BigNum>(basis)
}

/// Lattice reduction using the original Lenstra-Lenstra-Lovasz algorithm
///
/// This implementation uses platform double floating-point numbers (IEEE 754) for arithmetic operations.
/// The value of `delta` is set to 0.75.
///
///   - `basis`: A generating matrix for the lattice
///
/// The basis is reduced in-place.
#[deprecated(
    note = "Current implementation might yield incorrect results. Use l2.lll_float() instead"
)]
pub fn lll_float(basis: &mut Matrix<f64>) {
    lattice_reduce::<Float>(basis)
}
