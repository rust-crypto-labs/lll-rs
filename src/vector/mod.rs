//! Basic vector structures for LLL
use std::{
    fmt,
    ops::{self, Index, IndexMut},
};

pub type VectorF = Vector<f64>;
pub type BigVector = Vector<rug::Integer>;

/// Implementation of a vector without generic coefficients
#[derive(Clone)]
pub struct Vector<T> {
    /// Internal representation as a list of coefficients
    coefficients: Vec<T>,
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
    #![allow(dead_code)]
    fn basis_vector(dimension: usize, position: usize) -> Self {
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
    pub fn mulf(&self, other: &T) -> Self {
        let n = self.dimension();

        Self::from_vector(
            (0..n)
                .map(|i| self.coefficients[i].clone() * other)
                .collect(),
        )
    }
}

pub(crate) trait Dot {
    type Output;
    fn dot(&self, other: &Self) -> Self::Output;
}

impl Dot for BigVector {
    type Output = rug::Integer;
    fn dot(&self, other: &Self) -> Self::Output {
        self.coefficients
            .iter()
            .zip(&other.coefficients)
            .map(|(coeff_r, coeff_l)| coeff_r * coeff_l)
            .sum()
    }
}

impl Dot for VectorF {
    type Output = f64;
    fn dot(&self, other: &Self) -> Self::Output {
        self.coefficients
            .iter()
            .zip(&other.coefficients)
            .map(|(coeff_r, coeff_l)| coeff_r * coeff_l)
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
