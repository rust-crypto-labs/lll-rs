use crate::matrix::Matrix;
use crate::rug::{Integer, Rational};
use crate::vector::{BigVector, Dot, Vector, VectorF};
/**
 * Lattice reduction (Lenstra-Lenstra-Lovasz algorithm)
 */
pub fn lattice_reduce(basis: &mut Matrix<VectorF>) {
    // Parameter delta in the Lovasz condition
    let delta = 0.75;

    let (n, _) = basis.dimensions();
    let mut swap_condition = true;

    while swap_condition {
        // Perform rounded Gram-Schmidt orthogonalisation
        for i in 0..n {
            for k in 1..i {
                let j = i - k;

                let b_i = &basis[i];
                let b_j = &basis[j];
                let alpha = (b_i.dot(&b_j) / b_j.dot(&b_j)).round();
                basis[i] = b_i.sub(&b_j.mulf(alpha));
            }
        }

        // Check for the Lovasz condition and swap columns if appropriate
        swap_condition = false;
        for i in 0..n - 1 {
            let b_i = &basis[i];
            let b_ip1 = &basis[i + 1];

            let lhs = delta * b_i.dot(&b_i);

            let alpha = b_ip1.dot(&b_i) / b_i.dot(&b_i);
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

pub fn big_lattice_reduce(basis: &mut Matrix<BigVector>) {
    // Parameter delta in the Lovasz condition
    let delta = (3, 4);

    let (n, _) = basis.dimensions();
    let mut swap_condition = true;

    while swap_condition {
        // Perform rounded Gram-Schmidt orthogonalisation
        for i in 0..n {
            for k in 1..i {
                let j = i - k;

                let b_i = &basis[i];
                let b_j = &basis[j];
                let (_, alpha) =
                    Rational::from((b_i.dot(&b_j), b_j.dot(&b_j))).fract_round(Integer::new());
                basis[i] = b_i.sub(&b_j.mulf(alpha));
            }
        }

        // Check for the Lovasz condition and swap columns if appropriate
        swap_condition = false;
        for i in 0..n - 1 {
            let b_i = &basis[i];
            let b_ip1 = &basis[i + 1];

            let lhs = Rational::from(delta) * Rational::from(b_i.dot(&b_i));

            let (_, alpha) =
                Rational::from((b_ip1.dot(&b_i), b_i.dot(&b_i))).fract_round(Integer::new());
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
