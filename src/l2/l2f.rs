use crate::matrix::Matrix;
use crate::vector::{Dot, Vector, VectorF};

use std::cmp::max;

/**
 * Lattice reduction (L² algorithm)
 *
 * This implementation uses platform floating-point numbers (IEEE 754) for the
 * arithmetic operations.
 *
 * Arguments:
 *  * basis: A generating matrix for the lattice
 *  * eta:
 *  * delta:
 *
 * The basis is reduced in-place.
 * 
 * # Panics
 * if delta <= 1/4 or delta >= 1  
 * if eta <= 1/2 or eta > sqrt(delta)
 */
pub fn lattice_reduce(basis: &mut Matrix<VectorF>, eta: f64, delta: f64) {
    assert!(0.25<delta && delta<1.);
    assert!(0.5<eta && eta*eta<delta);
    // Variables
    let (d, _) = basis.dimensions();
    let mut gram: Matrix<VectorF> = Matrix::init(d, d); // Gram matrix (upper triangular)
    let mut r: Matrix<VectorF> = Matrix::init(d, d); // r_ij matrix
    let mut mu: Matrix<VectorF> = Matrix::init(d, d); // Gram coefficient matrix

    // Computing Gram matrix
    for i in 0..d {
        for j in 0..=i {
            gram[i][j] = basis[i].dot(&basis[j]);
        }
    }

    let eta_minus = (eta + 0.5) / 2.;
    let delta_plus = (delta + 1.) / 2.;

    r[0][0] = gram[0][0];

    let mut k = 1;

    while k < d {
        size_reduce(k, d, basis, &mut gram, &mut mu, &mut r, eta_minus);

        let delta_criterion = delta_plus * r[k - 1][k - 1];
        let scalar_criterion = r[k][k] + r[k - 1][k - 1] * mu[k][k - 1].powi(2);

        // Lovasz condition
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
                    r[i][j] =
                        gram[i][j] - (0..j).map(|index| mu[j][index] * r[i][index]).sum::<f64>();
                    mu[i][j] = r[i][j] / r[j][j];
                }
            }

            k = max(1, k - 1);
        }
    }
}

/**
 * Performs the `eta`-size-reduction of `basis[k]`
 *
 * Arguments:
 * * `k`: Index of the column to be `eta`-size-reduced
 * * `d`:
 * * `basis`: A generating
 * * `gram`: Gram matrix of `basis`  
 * * `mu`:
 * * `r`:
 * * `eta`:
 *
 * Note: both `basis` and `gram` are updated by this operation.
 */
fn size_reduce(
    k: usize,
    d: usize,
    basis: &mut Matrix<VectorF>,
    gram: &mut Matrix<VectorF>,
    mu: &mut Matrix<VectorF>,
    r: &mut Matrix<VectorF>,
    eta: f64,
) {
    // Update mu and r
    for i in 0..=k {
        r[k][i] = gram[k][i] - (0..i).map(|index| mu[i][index] * r[k][index]).sum::<f64>();
        mu[k][i] = r[k][i] / r[i][i];
    }

    if (0..k).any(|index| mu[k][index].abs() > eta) {
        for i in (0..k).rev() {
            let x = mu[k][i].round();
            basis[k] = basis[k].sub(&basis[i].mulf(x));

            // Updating Gram matrix
            for j in 0..d {
                if j < k {
                    gram[k][j] = basis[k].dot(&basis[j]);
                } else {
                    gram[j][k] = basis[k].dot(&basis[j]);
                }
            }

            for j in 0..i {
                mu[k][j] -= x * mu[i][j];
            }
        }
        size_reduce(k, d, basis, gram, mu, r, eta);
    }
}
