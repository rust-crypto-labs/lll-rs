//! Basic vector structures for LLL
use std::{
    fmt,
    ops::{self, Index, IndexMut},
};

pub type RationalVector = Vector<rug::Rational>;
pub type VectorF = Vector<f64>;
pub type BigVector = Vector<rug::Integer>;

/// Implementation of a vector with arbitrary-length rationals as
/// coefficients
#[derive(Clone)]
pub struct Vector<T> {
    /// Internal representation as a list of coefficients
    coefficients: Vec<T>,

    /// Dimension of the vector
    dimension: usize,
}

pub trait VectorMember:
    From<u32>
    + Clone
    + Default
    + for<'a> ops::Add<&'a Self, Output = Self>
    + for<'a> ops::Sub<&'a Self, Output = Self>
    + for<'a> std::ops::Mul<&'a Self, Output = Self>
    + std::iter::Sum<Self>
{
}

impl<T> VectorMember for T where
    T: From<u32>
        + Clone
        + Default
        + for<'a> ops::Add<&'a Self, Output = Self>
        + for<'a> ops::Sub<&'a Self, Output = Self>
        + for<'a> std::ops::Mul<&'a Self, Output = Self>
        + std::iter::Sum<Self>
{
}

impl<T> Vector<T>
where
    T: VectorMember,
{
    fn basis_vector(&self, position: usize) -> Self {
        assert!(position < self.dimension);

        let mut coefficients = vec![T::from(0); self.dimension()];
        coefficients[position] = T::from(1);

        Self {
            coefficients,
            dimension: self.dimension(),
        }
    }

    pub fn init(dimension: usize) -> Self {
        Self {
            coefficients: vec![Default::default(); dimension],
            dimension,
        }
    }

    pub fn dimension(&self) -> usize {
        self.dimension
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
        Self {
            dimension: coefficients.len(),
            coefficients,
        }
    }

    /// Multiplication by a scalar
    pub fn mulf(&self, other: &T) -> Self {
        let n = self.dimension();

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i].clone() * other)
                .collect(),
        )
    }

    pub fn dot(&self, other: &Self) -> T {
        let n = self.dimension();
        assert_eq!(n, other.dimension());

        (0..n)
            .map(|i| self.coefficients[i].clone() * &other.coefficients[i])
            .sum()
    }
}

impl<T> Index<usize> for Vector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &T {
        &self.coefficients[index]
    }
}

impl<T> IndexMut<usize> for Vector<T> {
    fn index_mut(&mut self, index: usize) -> &mut T {
        &mut self.coefficients[index]
    }
}

impl<T> fmt::Debug for Vector<T>
where
    T: fmt::Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self.coefficients)
    }
}
