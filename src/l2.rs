use crate::matrix::Matrix;
use crate::vector::{Dot, Vector, VectorF};

use std::cmp::max;

// TODO Algo in paper starts at 1 not 0 FIX IT

/**
 * Lattice reduction (LLL² algorithm)
 */
pub fn lattice_reduce(basis: &mut Matrix<VectorF>, mut eta: f64, mut delta: f64) {
    // Variables
    let d = basis.dimension;
    let mut gram: Matrix<VectorF> = Matrix::init(d); // Gram matix
    let mut r: Matrix<VectorF> = Matrix::init(d); // r_ij matrix
    let mut mu: Matrix<VectorF> = Matrix::init(d); // Gram coefficient matrix

    // Computing Gram matrix
    for i in 0..gram.dimension {
        for j in 0..gram.dimension {
            gram[i][j] = basis[i].dot(&basis[j]);
        }
    }

    eta += 0.5;
    eta /= 2.;
    delta += 1.;
    delta /= 2.;

    r[0][0] = gram[0][0];

    let mut k = 1;

    while k < d {
        size_reduce(k, eta, d, basis, &mut gram, &mut mu, &mut r);

        // Lovazs condition
        if delta * mu[k - 1][k - 1] < r[k][k] + mu[k][k - 1] * mu[k][k - 1] * r[k - 1][k - 1] {
            k += 1;
        } else {
            basis.swap(k, k - 1);

            // Updating Gram matrix
            for j in 0..d {
                gram[j][k] = basis[k].dot(&basis[j]);
                gram[k][j] = basis[k].dot(&basis[j]);
                gram[j][k - 1] = basis[k - 1].dot(&basis[j]);
                gram[k - 1][j] = basis[k - 1].dot(&basis[j]);
            }

            // Updating mu and r
            for i in 0..k {
                r[k][i] = gram[k][i] - (0..i).map(|index| mu[i][index] * r[k][index]).sum::<f64>();
                mu[k][i] = r[k][i] / r[i][i];
            }

            k = max(1, k - 1);
        }
    }
}

/**
 * Performs the eta-size-reduction of basis[k]
 */
fn size_reduce(
    k: usize,
    eta: f64,
    d: usize,
    basis: &mut Matrix<VectorF>,
    gram: &mut Matrix<VectorF>,
    mu: &mut Matrix<VectorF>,
    r: &mut Matrix<VectorF>,
) {
    // Update mu and r
    for i in 0..k {
        r[k][i] = gram[k][i] - (0..i).map(|index| mu[i][index] * r[k][index]).sum::<f64>();
        mu[k][i] = r[k][i] / r[i][i];
    }

    if (0..k).any(|index| mu[k][index].abs() > eta) {
        for i in (0..k).rev() {
            let x = mu[k][i].round();
            basis[k] = basis[k].sub(&basis[i].mulf(x));

            // Updating Gram matrix
            for j in 0..d {
                gram[j][k] = basis[k].dot(&basis[j]);
                gram[k][j] = basis[k].dot(&basis[j]);
            }

            for j in 0..i {
                mu[k][j] -= x * mu[i][j];
            }
            size_reduce(k, eta, d, basis, gram, mu, r);
        }
    }
}
