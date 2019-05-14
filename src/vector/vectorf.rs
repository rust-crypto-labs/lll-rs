use crate::vector::{Dot, Vector};

use std::{
    fmt,
    ops::{Index, IndexMut},
};

/**
 * Implementation of vectors with floating-point (IEEE 754) coefficients
 */
#[derive(Clone)]
pub struct VectorF {
    /// Underlying representation of the vector as a list of coefficients
    coefficients: Vec<f64>,

    /// Dimension of the vector
    dimension: usize,
}

impl Vector for VectorF {
    /**
     * Return a basis vector for the vector space
     *  `position`: number of the basis vector (0..n)
     */
    fn basis_vector(&self, position: usize) -> Self {
        assert!(position < self.dimension);

        let mut coefficients = vec![0.0; self.dimension()];
        coefficients[position] = 1.0;

        Self {
            coefficients,
            dimension: self.dimension(),
        }
    }

    /**
     * Create a new `VectorF` with default values, of size `dimension`
     */
    fn init(dimension: usize) -> Self {
        Self {
            coefficients: vec![Default::default(); dimension],
            dimension,
        }
    }

    /**
     * Return the vector's dimension
     */
    fn dimension(&self) -> usize {
        self.dimension
    }

    /**
     * Add two vectors of the same size
     */
    fn add(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i] + other.coefficients[i])
                .collect(),
        )
    }

    /**
     * Subtract the vector `other` from this vector
     */
    fn sub(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i] - other.coefficients[i])
                .collect(),
        )
    }
}

impl Dot<f64> for VectorF {
    /**
     * Dot product between two vectors
     */
    fn dot(&self, other: &Self) -> f64 {
        let n = self.dimension();
        assert_eq!(n, other.dimension());

        (0..n)
            .map(|i| self.coefficients[i] * other.coefficients[i])
            .sum()
    }
}

impl VectorF {
    /**
     * Create an instance from a `Vec` of floating-point coordinates
     */
    pub fn from_vector(coefficients: Vec<f64>) -> Self {
        Self {
            dimension: coefficients.len(),
            coefficients,
        }
    }

    /// Multiplication by a scalar
    pub fn mulf(&self, other: f64) -> Self {
        let n = self.dimension();

        Self::from_vector((0..n).map(|i| self.coefficients[i] * other).collect())
    }
}

impl Index<usize> for VectorF {
    type Output = f64;

    fn index(&self, index: usize) -> &f64 {
        &self.coefficients[index]
    }
}

impl IndexMut<usize> for VectorF {
    fn index_mut(&mut self, index: usize) -> &mut f64 {
        &mut self.coefficients[index]
    }
}

impl fmt::Debug for VectorF {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.coefficients)
    }
}
