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
//! let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
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

pub mod l2;
pub mod lll;
pub mod matrix;
pub mod vector;

#[cfg(test)]
mod test {
    use crate::{
        l2::{bigl2, l2f},
        lll::{biglll, lllf},
        matrix::Matrix,
        vector::{BigVector, VectorF},
    };

    use rug::{Assign, Integer};

    #[test]
    fn test_lllf() {
        // "Bad" lattice basis
        let mut basis: Matrix<VectorF> = Matrix::init(3, 4);
        basis[0] = VectorF::from_vector(vec![1., 0., 0., 1345.]);
        basis[1] = VectorF::from_vector(vec![0., 1., 0., 35.]);
        basis[2] = VectorF::from_vector(vec![0., 0., 1., 154.]);
        println!("{:?}", basis);

        // "Good" lattice basis
        lllf::lattice_reduce(&mut basis);
        println!("{:?}", basis);
    }

    #[test]
    fn test_biglll() {
        // "Bad" lattice basis
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0] = BigVector::from_vector(vec![
            Integer::from(1) << 100000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345),
        ]);
        basis[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(35),
        ]);
        basis[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(154),
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        biglll::lattice_reduce(&mut basis);
        println!("{:?}", basis);
    }

    #[test]
    fn test_l2f() {
        // "Bad" lattice basis
        let mut basis: Matrix<VectorF> = Matrix::init(3, 4);
        basis[0] = VectorF::from_vector(vec![1., 0., 0., 1345.]);
        basis[1] = VectorF::from_vector(vec![0., 1., 0., 35.]);
        basis[2] = VectorF::from_vector(vec![0., 0., 1., 154.]);
        println!("{:?}", basis);

        // "Good" lattice basis
        l2f::lattice_reduce(&mut basis, 0.501, 0.998);
        println!("{:?}", basis);
    }

    #[test]
    fn test_bigl2() {
        // "Bad" lattice basis
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0] = BigVector::from_vector(vec![
            Integer::from(1) << 100000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345),
        ]);
        basis[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(35),
        ]);
        basis[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(154),
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        bigl2::lattice_reduce(&mut basis, 0.5005, 0.999);
        println!("{:?}", basis);
    }

    #[test]
    fn test_bigl2_2() {
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0][0].assign(1);
        basis[0][3].assign(1);
        basis[1][1].assign(1);
        basis[1][3].assign(5);
        basis[2][2].assign(1);
        basis[2][3].assign(9);
        bigl2::lattice_reduce(&mut basis, 0.6, 0.95);
    }
}
