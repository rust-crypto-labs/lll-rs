use crate::matrix::Matrix;
use crate::scalars::{Scalars, FromExt};
use crate::vector::{Dot, Vector, Coefficient};

use std::cmp::max;

/// Lattice reduction (L² algorithm)
///
/// This implementation uses generic Scalar types for the underlying arithmetic operations.
///
/// Arguments:
///  * basis: A generating matrix for the lattice
///  * eta: eta factor of the basis reduction
///  * delta: delta factor of the basis reduction
///
/// The basis is reduced in-place.
///
/// # Panics
/// if delta <= 1/4 or delta >= 1  
/// if eta <= 1/2 or eta > sqrt(delta)
pub(crate) fn lattice_reduce<S>(basis: &mut Matrix<S::Integer>, eta: f64, delta: f64)
where
    S: Scalars,
    S::Integer: Coefficient,
    S::Fraction: Coefficient,
    Vector<S::Integer>: Dot<Output = S::Integer>,
{
    assert!(0.25 < delta && delta < 1.);
    assert!(0.5 < eta && eta * eta < delta);
    // Variables
    let (d, _) = basis.dimensions();
    let mut gram: Matrix<S::Integer> = Matrix::init(d, d); // Gram matrix (upper triangular)
    let mut r: Matrix<S::Fraction> = Matrix::init(d, d); // r_ij matrix
    let mut mu: Matrix<S::Fraction> = Matrix::init(d, d); // Gram coefficient matrix

    // Computing Gram matrix
    for i in 0..d {
        for j in 0..=i {
            gram[i][j] = basis[i].dot(&basis[j]);
        }
    }

    let eta_minus = S::Fraction::from_ext((eta + 0.5) / 2.);
    let delta_plus = S::Fraction::from_ext((delta + 1.) / 2.);

    r[0][0] = S::Fraction::from_ext(&gram[0][0]);

    let mut k = 1;

    while k < d {
        size_reduce::<S>(k, d, basis, &mut gram, &mut mu, &mut r, &eta_minus);

        let delta_criterion = delta_plus.clone() * &r[k - 1][k - 1];
        let scalar_criterion =
            (mu[k][k - 1].clone() * &mu[k][k - 1] * &r[k - 1][k - 1]) + &r[k][k];

        // Lovazs condition
        if delta_criterion < scalar_criterion {
            k += 1;
        } else {
            basis.swap(k, k - 1);

            // Updating Gram matrix
            for j in 0..d {
                if j < k {
                    gram[k][j] = basis[k].dot(&basis[j]);
                    gram[k - 1][j] = basis[k - 1].dot(&basis[j]);
                } else {
                    gram[j][k] = basis[k].dot(&basis[j]);
                    gram[j][k - 1] = basis[k - 1].dot(&basis[j]);
                }
            }

            // Updating mu and r
            for i in 0..=k {
                for j in 0..=i {
                    r[i][j] = S::Fraction::from_ext(&gram[i][j])
                        - &(0..j)
                            .map(|index| mu[j][index].clone() * &r[i][index])
                            .sum::<S::Fraction>();
                    mu[i][j] = r[i][j].clone() / &r[j][j];
                }
            }

            k = max(1, k - 1);
        }
    }
}

/// Performs the `eta`-size-reduction of `basis[k]`
///
/// Arguments:
/// * `k`: Index of the column to be `eta`-size-reduced
/// * `d`: The basis dimension
/// * `basis`: A generating matrix for the lattice
/// * `gram`: Gram matrix of `basis`  
/// * `mu`: Gram coefficient matrix
/// * `r`: the r_ij matrix
/// * `eta`: eta factor of the basis reduction
///
/// Note: both `basis` and `gram` are updated by this operation.
fn size_reduce<S>(
    k: usize,
    d: usize,
    basis: &mut Matrix<S::Integer>,
    gram: &mut Matrix<S::Integer>,
    mu: &mut Matrix<S::Fraction>,
    r: &mut Matrix<S::Fraction>,
    eta: &S::Fraction,
) where
    S: Scalars,
    S::Integer: Coefficient,
    S::Fraction: Coefficient,
    Vector<S::Integer>: Dot<Output = S::Integer>,
{
    // Update mu and r
    for i in 0..=k {
        r[k][i] = S::Fraction::from_ext(&gram[k][i])
            - &(0..i)
                .map(|index| mu[i][index].clone() * &r[k][index])
                .sum::<S::Fraction>();
        mu[k][i] = r[k][i].clone() / &r[i][i];
    }

    if (0..k).any(|index| S::abs(mu[k][index].clone()) > *eta) {
        for i in (0..k).rev() {
            let x = S::round(&mu[k][i]);
            basis[k] = basis[k].sub(&basis[i].mulf(&x));

            // Updating Gram matrix
            for j in 0..d {
                if j < k {
                    gram[k][j] = basis[k].dot(&basis[j]);
                } else {
                    gram[j][k] = basis[k].dot(&basis[j]);
                }
            }

            for j in 0..i {
                let minus = S::Fraction::from_ext(&x) * &mu[i][j];
                mu[k][j] -= &minus;
            }
        }
        size_reduce::<S>(k, d, basis, gram, mu, r, eta);
    }
}

pub mod bigl2 {
    use crate::matrix::Matrix;
    use crate::scalars::BigNum;

    /// Lattice reduction (L² algorithm)
    ///
    /// This implementation uses `rug::Integers` and `rug::Rationnal` for the underlying arithmetic operations.
    ///
    /// Arguments:
    ///  * basis: A generating matrix for the lattice
    ///  * eta: eta factor of the basis reduction
    ///  * delta: delta factor of the basis reduction
    ///
    /// The basis is reduced in-place.
    ///
    /// # Panics
    /// if delta <= 1/4 or delta >= 1  
    /// if eta <= 1/2 or eta > sqrt(delta)
    pub fn lattice_reduce(basis: &mut Matrix<rug::Integer>, eta: f64, delta: f64) {
        super::lattice_reduce::<BigNum>(basis, eta, delta)
    }
}

pub mod l2f {
    use crate::matrix::Matrix;
    use crate::scalars::Float;

    /// Lattice reduction (L² algorithm)
    ///
    /// This implementation uses platform double floating-point numbers (IEEE 754)
    /// for the underlying arithmetic operations.
    ///
    /// Arguments:
    ///  * basis: A generating matrix for the lattice
    ///  * eta: eta factor of the basis reduction
    ///  * delta: delta factor of the basis reduction
    ///
    /// The basis is reduced in-place.
    ///
    /// # Panics
    /// if delta <= 1/4 or delta >= 1  
    /// if eta <= 1/2 or eta > sqrt(delta)
    pub fn lattice_reduce(basis: &mut Matrix<f64>, eta: f64, delta: f64) {
        super::lattice_reduce::<Float>(basis, eta, delta)
    }
}
