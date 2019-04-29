use fmt::Debug;
use std::fmt;
use std::ops::{Index, IndexMut};

use crate::vector::Vector;

/// Matrix
pub struct Matrix<T: Vector> {
    columns: Vec<T>,
    dimensions: (usize, usize),
}

impl<T: Vector> Matrix<T>
where
    T: Clone,
{
    /// Initialise matrix type
    pub fn init(col_num: usize, col_dim: usize) -> Self {
        Self {
            columns: vec![T::init(col_dim); col_num],
            dimensions: (col_num, col_dim),
        }
    }

    /**
     * Return the matrix dimensions
     */
    pub fn dimensions(&self) -> (usize, usize) {
        self.dimensions
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
