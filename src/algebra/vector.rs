//! Basic vector structures for LLL
use super::Coefficient;

use std::{
    fmt,
    ops::{Index, IndexMut},
};

/// Implementation of a vector without generic coefficients
#[derive(Clone, PartialEq)]
pub struct Vector<T: Coefficient> {
    /// Internal representation as a list of coefficients
    coefficients: Vec<T>,
}

impl<T: Coefficient> Vector<T> {
    pub fn basis_vector(dimension: usize, position: usize) -> Self {
        assert!(position < dimension);

        let coefficients = (0..dimension)
            .map(|i| {
                if i == position {
                    T::from(1)
                } else {
                    T::from(0)
                }
            })
            .collect();

        Self { coefficients }
    }

    pub fn init(dimension: usize) -> Self {
        Self {
            coefficients: vec![Default::default(); dimension],
        }
    }

    pub fn dimension(&self) -> usize {
        self.coefficients.len()
    }

    pub fn add(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i].clone() + &other.coefficients[i])
                .collect(),
        )
    }

    pub fn sub(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i].clone() - &other.coefficients[i])
                .collect(),
        )
    }

    /// Create an instance from a `Vec`
    pub fn from_vector(coefficients: Vec<T>) -> Self {
        Self { coefficients }
    }

    /// Multiplication by a scalar
    pub fn mulf(&self, other: T) -> Self {
        let n = self.dimension();

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i].clone() * &other)
                .collect(),
        )
    }
    pub fn zero(dimension: usize) -> Self {
        Self {
            coefficients: vec![Default::default(); dimension],
        }
    }

    pub fn is_zero(&self) -> bool {
        self == &Vector::zero(self.dimension())
    }

    pub fn dot(&self, other: &Self) -> T {
        self.coefficients
            .iter()
            .zip(&other.coefficients)
            .map(|(coeff_r, coeff_l)| coeff_r.clone() * coeff_l)
            .sum()
    }
}

impl<T: Coefficient> Index<usize> for Vector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &T {
        &self.coefficients[index]
    }
}

impl<T: Coefficient> IndexMut<usize> for Vector<T> {
    fn index_mut(&mut self, index: usize) -> &mut T {
        &mut self.coefficients[index]
    }
}

impl<T: Coefficient> fmt::Debug for Vector<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{:?}", self.coefficients)
    }
}
/*
pub type VectorF = Vector<Float>;
pub type BigVector = Vector<BigNum>;*/
