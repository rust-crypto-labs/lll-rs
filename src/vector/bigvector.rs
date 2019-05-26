use crate::vector::{Dot, Vector};

use rug::Integer;

use std::{
    fmt,
    ops::{Index, IndexMut},
};

/**
 * Implementation of vectors with arbitrary-length integers as underlying
 * coefficients
 */
#[derive(Clone)]
pub struct BigVector {
    /// Internal representation as a list of `Integer`s
    coefficients: Vec<Integer>,

    /// Dimension of the vector
    dimension: usize,
}

impl Vector for BigVector {
    fn basis_vector(&self, position: usize) -> Self {
        assert!(position < self.dimension);

        let mut coefficients = vec![Integer::from(0); self.dimension()];
        coefficients[position] = Integer::from(1);

        Self {
            coefficients,
            dimension: self.dimension(),
        }
    }

    fn init(dimension: usize) -> Self {
        Self {
            coefficients: vec![Default::default(); dimension],
            dimension,
        }
    }

    fn dimension(&self) -> usize {
        self.dimension
    }

    fn add(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| Integer::from(&self.coefficients[i] + &other.coefficients[i]))
                .collect(),
        )
    }

    fn sub(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| Integer::from(&self.coefficients[i] - &other.coefficients[i]))
                .collect(),
        )
    }
}

impl BigVector {
    /**
     * Create an instance from a `Vec` of `Integer`s
     */
    pub fn from_vector(coefficients: Vec<Integer>) -> Self {
        Self {
            dimension: coefficients.len(),
            coefficients,
        }
    }

    /// Multiplication by a scalar
    pub fn mulf(&self, other: &Integer) -> Self {
        let n = self.dimension();

        Self::from_vector(
            (0..n)
                .map(|i| Integer::from(&self.coefficients[i] * other))
                .collect(),
        )
    }
}

impl Dot<Integer> for BigVector {
    fn dot(&self, other: &Self) -> Integer {
        let n = self.dimension();
        assert_eq!(n, other.dimension());

        (0..n)
            .map(|i| Integer::from(&self.coefficients[i] * &other.coefficients[i]))
            .sum()
    }
}

impl Index<usize> for BigVector {
    type Output = Integer;

    fn index(&self, index: usize) -> &Integer {
        &self.coefficients[index]
    }
}

impl IndexMut<usize> for BigVector {
    fn index_mut(&mut self, index: usize) -> &mut Integer {
        &mut self.coefficients[index]
    }
}

impl fmt::Debug for BigVector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.coefficients)
    }
}
