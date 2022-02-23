use crate::algebra::{BigNum, Float, FromExt, Matrix, Scalar, Vector};

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
fn lattice_reduce<S: Scalar>(basis: &mut Matrix<S::Integer>, eta: f64, delta: f64) {
    assert!(0.25 < delta && delta < 1.);
    assert!(0.5 < eta && eta * eta < delta);

    // Variables
    let (d, _) = basis.dimensions();
    let mut gram: Matrix<S::Integer> = Matrix::init(d, d); // Gram matrix (upper triangular)
    let mut r: Matrix<S::Fraction> = Matrix::init(d, d); // r_ij matrix
    let mut mu: Matrix<S::Fraction> = Matrix::init(d, d); // Gram coefficient matrix
    let mut s: Vector<S::Fraction> = Vector::init(d);

    let zero = S::Fraction::from(0);
    let mut zeros = 0;

    // Computing Gram matrix
    for i in 0..d {
        for j in 0..=i {
            gram[i][j] = basis[i].dot(&basis[j]);
        }
    }

    let eta_minus = S::Fraction::from_ext((eta + 0.5) / 2.);
    let delta_plus = S::Fraction::from_ext(0.99); //(delta + 1.) / 2.);

    r[0][0] = S::Fraction::from_ext(&gram[0][0]);

    let mut kappa = 1;

    while kappa < (d - zeros) {
        println!("Before size reduce:");
        println!("{:?}", &basis);
        println!("-------------------");
        size_reduce::<S>(basis, &mut gram, &mut mu, &mut r, kappa, &eta_minus);

        s[0] = S::Fraction::from_ext((gram[kappa][kappa].clone(), S::Integer::from(1)));
        for i in 0..kappa {
            s[i + 1] = s[i].clone() - &(mu[kappa][i].clone() * &r[kappa][i]);
        }

        let delta_criterion = delta_plus.clone() * &r[kappa - 1][kappa - 1];

        if delta_criterion > s[kappa - 1] {
            let kappa_prime = kappa;
            while kappa >= 1 && delta_criterion >= s[kappa - 1] {
                kappa -= 1;
            }

            let is_neg = s[kappa] <= zero;

            let k = if !is_neg {
                kappa
            } else {
                zeros += 1;
                d - zeros
            };

            if k != kappa_prime {
                basis.insert(kappa_prime, k);
                mu.insert(kappa_prime, k);
                r.insert(kappa_prime, k)
            }

            // Update Gram matrix
            for i in 0..d {
                for j in 0..=i {
                    gram[i][j] = basis[i].dot(&basis[j]);
                }
            }

            if is_neg {
                kappa = kappa_prime;
                continue;
            }
        }
        r[kappa][kappa] = s[kappa].clone();
        kappa += 1;
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
fn size_reduce<S: Scalar>(
    basis: &mut Matrix<S::Integer>,
    gram: &mut Matrix<S::Integer>,
    mu: &mut Matrix<S::Fraction>,
    r: &mut Matrix<S::Fraction>,
    kappa: usize,
    eta: &S::Fraction,
) {
    let zero = S::Integer::from(0);
    let one = S::Integer::from(1);
    loop {
        cfa::<S>(kappa, basis, gram, mu, r);

        let all_zeroes = (0..kappa)
            .rev()
            .all(|i| &S::abs(mu[kappa][i].clone()) < eta);

        if all_zeroes {
            break;
        }

        let mut m = mu[kappa].clone();

        for i in (0..kappa).rev() {
            let x_i = S::round(&m[i]);
            if x_i != zero {
                for j in 0..i {
                    m[j] -=
                        &(mu[i][j].clone() * &S::Fraction::from_ext((x_i.clone(), one.clone())));
                }

                // Swap basis
                basis[kappa] = basis[kappa].sub(&basis[i].mulf(x_i));
            }
        }

        // Update Gram matrix

        for j in 0..=kappa {
            gram[kappa][j] = basis[kappa].dot(&basis[j]);
        }
    }
}

fn cfa<S: Scalar>(
    i: usize,
    basis: &mut Matrix<S::Integer>,
    gram: &mut Matrix<S::Integer>,
    mu: &mut Matrix<S::Fraction>,
    r: &mut Matrix<S::Fraction>,
) {
    for j in 0..=i {
        gram[i][j] = basis[i].dot(&basis[j]);
    }

    for j in 0..i {
        r[i][j] = S::Fraction::from_ext((gram[i][j].clone(), S::Integer::from(1)));

        for k in 0..j {
            r[i][j] = r[i][j].clone() - &(r[i][k].clone() * &mu[j][k]);
        }
        mu[i][j] = r[i][j].clone() / &r[j][j];
    }
}

/// Puts the trailing null columns at the beginning of the matrix
fn zeros_first<S: Scalar>(basis: &mut Matrix<S::Integer>) {
    let (d, _) = basis.dimensions();
    while basis[d - 1].is_zero() {
        basis.insert(d - 1, 0)
    }
}

fn reduction<S: Scalar>(basis: &mut Matrix<S::Integer>, eta: f64, delta: f64) {
    lattice_reduce::<S>(basis, eta, delta);
    lattice_reduce::<S>(basis, eta, delta);
    zeros_first::<S>(basis);
}

/// Lattice reduction (L² algorithm)
///
/// This implementation uses `rug::Integers` and `rug::Rationnal` for the underlying arithmetic operations.
///
/// Arguments:
///  * basis: A generating matrix for the lattice
///  * eta: eta factor of the basis reduction
///  * delta: delta factor of the basis reduction
///
/// The basis is reduced in-place. The reduction is performed according to the standard pipeline of the fplll implementation of LLL.
/// It is done by doing one extra LLL-reduction at the end and putting all the trailing null rows at the beginning
///
/// # Panics
/// if delta <= 1/4 or delta >= 1  
/// if eta <= 1/2 or eta > sqrt(delta)
pub fn lll_bignum(basis: &mut Matrix<rug::Integer>, eta: f64, delta: f64) {
    reduction::<BigNum>(basis, eta, delta)
}

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
/// The basis is reduced in-place. The reduction is performed according to the standard pipeline of the fplll implementation of LLL.
/// It is done by doing one extra LLL-reduction at the end and putting all the trailing null rows at the beginning
///
/// # Panics
/// if delta <= 1/4 or delta >= 1  
/// if eta <= 1/2 or eta > sqrt(delta)
pub fn lll_float(basis: &mut Matrix<f64>, eta: f64, delta: f64) {
    reduction::<Float>(basis, eta, delta)
}
