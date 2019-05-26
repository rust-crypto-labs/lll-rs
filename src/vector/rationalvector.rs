use crate::vector::{Dot, Vector};

use rug::Rational;

use std::{
    fmt,
    ops::{Index, IndexMut},
};

/**
 * Implementation of a vector with arbitrary-length rationals as
 * coefficients
 */
#[derive(Clone)]
pub struct RationalVector {
    /// Internal representation as a list of coefficients
    coefficients: Vec<Rational>,

    /// Dimension of the vector
    dimension: usize,
}

impl Vector for RationalVector {
    fn basis_vector(&self, position: usize) -> Self {
        assert!(position < self.dimension);

        let mut coefficients = vec![Rational::from(0); self.dimension()];
        coefficients[position] = Rational::from(1);

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
                .map(|i| Rational::from(&self.coefficients[i] + &other.coefficients[i]))
                .collect(),
        )
    }

    fn sub(&self, other: &Self) -> Self {
        let n = self.dimension();

        assert_eq!(n, other.dimension());

        Self::from_vector(
            (0..n)
                .map(|i| Rational::from(&self.coefficients[i] - &other.coefficients[i]))
                .collect(),
        )
    }
}

impl RationalVector {
    /**
     * Create an instance from a `Vec` of `Rational`
     */
    pub fn from_vector(coefficients: Vec<Rational>) -> Self {
        Self {
            dimension: coefficients.len(),
            coefficients,
        }
    }

    /// Multiplication by a scalar
    pub fn mulf(&self, other: &Rational) -> Self {
        let n = self.dimension();

        Self::from_vector(
            (0..n)
                .map(|i| Rational::from(&self.coefficients[i] * other))
                .collect(),
        )
    }
}

impl Dot<Rational> for RationalVector {
    fn dot(&self, other: &Self) -> Rational {
        let n = self.dimension();
        assert_eq!(n, other.dimension());

        (0..n)
            .map(|i| Rational::from(&self.coefficients[i] * &other.coefficients[i]))
            .sum()
    }
}

impl Index<usize> for RationalVector {
    type Output = Rational;

    fn index(&self, index: usize) -> &Rational {
        &self.coefficients[index]
    }
}

impl IndexMut<usize> for RationalVector {
    fn index_mut(&mut self, index: usize) -> &mut Rational {
        &mut self.coefficients[index]
    }
}

impl fmt::Debug for RationalVector {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.coefficients)
    }
}
