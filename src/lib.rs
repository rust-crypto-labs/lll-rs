//! A Rust implementation of the (basic) Lenstra-Lenstra-Lovasz lattice reduction algorithm
//!
//! # Introduction
//! `lll-rs` is an implementation of the Lenstra-Lenstra-Lovász lattice basis reduction
//! algorithm in Rust as well as the implementation of an improved version, the L² algorithm.
//! The library comes with a set of simple helpers to create vectors and matrices to perform
//! lattice basis reduction.
//!
//! # Examples
//!
//! ```rust
//! use lll_rs::{
//!     l2::{bigl2, l2f},
//!     lll::{biglll, lllf},
//!     matrix::Matrix,
//!     vector::{BigVector, VectorF},
//! };
//!
//! use rug::{Integer,Assign};
//!
//! // Init the matrix with Integer
//! let mut basis: Matrix<Integer> = Matrix::init(3, 4);
//!
//! // Populate the matix
//! basis[0] = BigVector::from_vector(vec![
//!     Integer::from(1) << 100000,
//!     Integer::from(0),
//!     Integer::from(0),
//!     Integer::from(1345),
//! ]);
//! basis[1] = BigVector::from_vector(vec![
//!     Integer::from(0),
//!     Integer::from(1),
//!     Integer::from(0),
//!     Integer::from(35),
//! ]);
//! basis[2] = BigVector::from_vector(vec![
//!     Integer::from(0),
//!     Integer::from(0),
//!     Integer::from(1),
//!     Integer::from(154),
//! ]);
//!
//! // Perfom the LLL basis redution
//! biglll::lattice_reduce(&mut basis);
//!
//! // OR
//! // Perfom the LLL basis redution
//! // Specify the delta and eta coefficient for the reduction
//! bigl2::lattice_reduce(&mut basis, 0.5005, 0.999);
//! ```
//!
extern crate rug;

mod algebra;
pub mod l2;
pub mod lll;

pub use algebra::{BigNum, Float, Matrix};

#[cfg(test)]
mod test {
    use crate::{l2, lll, Matrix};

    #[test]
    fn test_lllf() {
        // "Bad" lattice basis
        let mut basis: Matrix<f64> = Matrix::from_matrix(vec![
            vec![1., 0., 0., 1345.],
            vec![0., 1., 0., 35.],
            vec![0., 0., 1., 154.],
        ]);

        // "Good" lattice basis
        lll::lll_float(&mut basis);

        let result: Matrix<f64> = Matrix::from_matrix(vec![
            vec![0.0, -4.0, 1.0, 14.0],
            vec![0.0, 1.0, 0.0, 35.0],
            vec![1.0, 348.0, -88.0, -27.0],
        ]);

        assert_eq!(basis, result);
    }

    #[test]
    fn test_biglll() {
        type I = rug::Integer;
        // "Bad" lattice basis
        let mut basis: Matrix<I> = Matrix::from_matrix(vec![
            vec![I::from(1) << 100000, I::from(0), I::from(0), I::from(1345)],
            vec![I::from(0), I::from(1), I::from(0), I::from(35)],
            vec![I::from(0), I::from(0), I::from(1), I::from(154)],
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        lll::lll_bignum(&mut basis);
        println!("{:?}", basis);
    }

    #[test]
    fn test_l2f() {
        // "Bad" lattice basis
        let mut basis: Matrix<f64> = Matrix::from_matrix(vec![
            vec![1., 0., 0., 1345.],
            vec![0., 1., 0., 35.],
            vec![0., 0., 1., 154.],
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        l2::lll_float(&mut basis, 0.501, 0.998);
        println!("{:?}", basis);

        let result: Matrix<f64> = Matrix::from_matrix(vec![
            vec![1.0, 1.0, -9.0, -6.0],
            vec![0.0, 9.0, -2.0, 7.0],
            vec![1.0, -3.0, -8.0, 8.0],
        ]);

        assert_eq!(basis, result);
    }

    #[test]
    fn test_bigl2() {
        type I = rug::Integer;
        let mut basis: Matrix<I> = Matrix::from_matrix(vec![
            vec![I::from(1), I::from(2), I::from(3)],
            vec![I::from(4), I::from(5), I::from(6)],
            vec![I::from(7), I::from(8), I::from(9)],
        ]);
        println!("{:?}", basis);

        l2::lll_bignum(&mut basis, 0.6, 0.95);
        println!("{:?}", basis);

        let result: Matrix<I> = Matrix::from_matrix(vec![
            vec![I::from(0), I::from(0), I::from(0)],
            vec![I::from(2), I::from(1), I::from(0)],
            vec![I::from(-1), I::from(1), I::from(3)],
        ]);
        assert_eq!(basis, result);
    }

    #[test]
    fn test_bigl2_ntrulike() {
        type I = rug::Integer;
        let mut basis: Matrix<I> = Matrix::from_matrix(vec![
            vec![
                I::from(1),
                I::from(0),
                I::from(0),
                I::from(436),
                I::from(225),
                I::from(381),
            ],
            vec![
                I::from(0),
                I::from(1),
                I::from(0),
                I::from(381),
                I::from(436),
                I::from(225),
            ],
            vec![
                I::from(0),
                I::from(0),
                I::from(1),
                I::from(225),
                I::from(381),
                I::from(436),
            ],
            vec![
                I::from(0),
                I::from(0),
                I::from(0),
                I::from(521),
                I::from(0),
                I::from(0),
            ],
            vec![
                I::from(0),
                I::from(0),
                I::from(0),
                I::from(0),
                I::from(521),
                I::from(0),
            ],
            vec![
                I::from(0),
                I::from(0),
                I::from(0),
                I::from(0),
                I::from(0),
                I::from(521),
            ],
        ]);
        println!("{:?}", basis);

        l2::lll_bignum(&mut basis, 0.6, 0.95);
        println!("{:?}", basis);

        let result: Matrix<I> = Matrix::from_matrix(vec![
            vec![
                I::from(1),
                I::from(1),
                I::from(1),
                I::from(0),
                I::from(0),
                I::from(0),
            ],
            vec![
                I::from(-11),
                I::from(0),
                I::from(12),
                I::from(-12),
                I::from(13),
                I::from(-1),
            ],
            vec![
                I::from(12),
                I::from(-11),
                I::from(0),
                I::from(-1),
                I::from(-12),
                I::from(13),
            ],
            vec![
                I::from(12),
                I::from(-1),
                I::from(-10),
                I::from(-4),
                I::from(17),
                I::from(-13),
            ],
            vec![
                I::from(1),
                I::from(10),
                I::from(-12),
                I::from(-17),
                I::from(13),
                I::from(4),
            ],
            vec![
                I::from(8),
                I::from(-5),
                I::from(-4),
                I::from(162),
                I::from(180),
                I::from(179),
            ],
        ]);
        assert_eq!(basis, result);
    }
}
