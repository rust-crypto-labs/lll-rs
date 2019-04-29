use crate::matrix::Matrix;
use crate::rug::{Integer, Rational};
use crate::vector::{BigVector, Dot, RationalVector, Vector, VectorF};

use std::cmp::max;
use std::cmp::Ordering;

/**
 * Lattice reduction (LLL² algorithm)
 */
pub fn lattice_reduce(basis: &mut Matrix<VectorF>, eta: f64, delta: f64) {
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

        // Lovazs condition
        if delta_plus * r[k - 1][k - 1] < r[k][k] + mu[k][k - 1] * mu[k][k - 1] * r[k - 1][k - 1] {
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
 * Performs the eta-size-reduction of basis[k]
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

pub fn big_lattice_reduce(basis: &mut Matrix<BigVector>, eta: f64, delta: f64) {
    // Variables
    let (d, _) = basis.dimensions();
    let mut gram: Matrix<BigVector> = Matrix::init(d, d); // Gram matrix (upper triangular)
    let mut r: Matrix<RationalVector> = Matrix::init(d, d); // r_ij matrix
    let mut mu: Matrix<RationalVector> = Matrix::init(d, d); // Gram coefficient matrix

    // Computing Gram matrix
    for i in 0..d {
        for j in 0..=i {
            gram[i][j] = basis[i].dot(&basis[j]);
        }
    }

    let eta_minus = Rational::from_f64((eta + 0.5) / 2.).unwrap();
    let delta_plus = Rational::from_f64((delta + 1.) / 2.).unwrap();

    r[0][0] = Rational::from(&gram[0][0]);

    let mut k = 1;

    while k < d {
        big_size_reduce(k, d, basis, &mut gram, &mut mu, &mut r, &eta_minus);

        let delta_criterion = Rational::from(&delta_plus * &r[k - 1][k - 1]);
        let scalar_criterion = &r[k][k] + Rational::from(&mu[k][k - 1] * &r[k - 1][k - 1]);

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
                    r[i][j] = Rational::from(&gram[i][j])
                        - (0..j)
                            .map(|index| Rational::from(&mu[j][index]).square() * &r[i][index])
                            .sum::<Rational>();
                    mu[i][j] = Rational::from(&r[i][j] / &r[j][j]);
                }
            }

            k = max(1, k - 1);
        }
    }
}

/**
 * Performs the eta-size-reduction of basis[k]
 */
fn big_size_reduce(
    k: usize,
    d: usize,
    basis: &mut Matrix<BigVector>,
    gram: &mut Matrix<BigVector>,
    mu: &mut Matrix<RationalVector>,
    r: &mut Matrix<RationalVector>,
    eta: &Rational,
) {
    // Update mu and r
    for i in 0..=k {
        r[k][i] = Rational::from(&gram[k][i])
            - (0..i)
                .map(|index| Rational::from(&mu[i][index] * &r[k][index]))
                .sum::<Rational>();
        mu[k][i] = Rational::from(&r[k][i] / &r[i][i]);
    }

    if (0..k).any(|index| eta.cmp_abs(&mu[k][index]) == Ordering::Less) {
        for i in (0..k).rev() {
            let (_, x) = mu[k][i].clone().fract_round(Integer::new());
            basis[k] = basis[k].sub(&basis[i].mulf(x.clone()));

            // Updating Gram matrix
            for j in 0..d {
                if j < k {
                    gram[k][j] = basis[k].dot(&basis[j]);
                } else {
                    gram[j][k] = basis[k].dot(&basis[j]);
                }
            }

            for j in 0..i {
                let shift = mu[i][j].clone();
                mu[k][j] -= Rational::from(&x) * shift;
            }
        }
        big_size_reduce(k, d, basis, gram, mu, r, eta);
    }
}
