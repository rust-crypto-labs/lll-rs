use fmt::Debug;
use std::fmt;
use std::ops::{Index, IndexMut};

use crate::vector::Vector;

/// Matrix
pub struct Matrix<T: Vector> {
    pub columns: Vec<T>,
    pub dimension: usize,
}

impl<T: Vector> Matrix<T> {
    /// Initialise matrix type
    pub fn init(dimension: usize) -> Self {
        Self {
            columns: vec![],
            dimension,
        }
    }

    /// Identity matrix
    pub fn identity(&self) -> Self {
        let n = self.dimension;
        let vector_type = T::init(n);

        Self {
            columns: (0..n).map(|i| vector_type.basis_vector(i)).collect(),
            dimension: n,
        }
    }

    pub fn swap(&mut self, i: usize, j: usize) {
        self.columns.swap(i, j);
    }
}

impl<T> Index<usize> for Matrix<T>
where
    T: Vector,
{
    type Output = T;

    fn index(&self, index: usize) -> &T {
        &self.columns[index]
    }
}

impl<T> IndexMut<usize> for Matrix<T>
where
    T: Vector,
{
    fn index_mut(&mut self, index: usize) -> &mut T {
        &mut self.columns[index]
    }
}

impl<T> fmt::Debug for Matrix<T>
where
    T: Vector + Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{:?}\n", self.columns)
    }
}
