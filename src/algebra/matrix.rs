//! Basic matrix structure for LLL

use super::{Coefficient, Vector};

use std::{
    fmt,
    ops::{Index, IndexMut},
};

#[derive(PartialEq)]
/// A `Matrix` is a collection of `Vector`s
pub struct Matrix<T: Coefficient> {
    /// Internal representation as a list of elements of type `T`
    columns: Vec<Vector<T>>,

    /// Dimensions of the matrix
    dimensions: (usize, usize),
}

impl<T: Coefficient> Matrix<T> {
    /// Initialise an empty `Matrix`
    ///      - `col_num`: number of columns
    ///      - `col_dim`: number of rows
    pub fn init(col_num: usize, col_dim: usize) -> Self {
        Self {
            columns: vec![Vector::<T>::init(col_dim); col_num],
            dimensions: (col_num, col_dim),
        }
    }

    pub fn from_columns(columns: Vec<Vector<T>>) -> Self {
        let dimensions = if let Some(col) = columns.first() {
            (columns.len(), col.dimension())
        } else {
            (0, 0)
        };
        Self {
            columns,
            dimensions,
        }
    }

    pub fn from_matrix(matrix: Vec<Vec<T>>) -> Self {
        Self::from_columns(
            matrix
                .iter()
                .map(|column| Vector::<T>::from_vector(column.to_vec()))
                .collect(),
        )
    }

    /// Return the matrix dimensions
    pub fn dimensions(&self) -> (usize, usize) {
        self.dimensions
    }

    /// Swap two columns of the matrix
    pub fn swap(&mut self, i: usize, j: usize) {
        self.columns.swap(i, j);
    }

    /// Insert the i-th column before the j-th one
    pub fn insert(&mut self, i: usize, j: usize) {
        let v = self.columns.remove(i);
        self.columns.insert(j, v)
    }
}

/// Direct access to a column
impl<T: Coefficient> Index<usize> for Matrix<T> {
    type Output = Vector<T>;

    fn index(&self, index: usize) -> &Self::Output {
        &self.columns[index]
    }
}

/// Direct access to a column (mutable)
impl<T: Coefficient> IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.columns[index]
    }
}

impl<T: Coefficient> fmt::Debug for Matrix<T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "{:?}", self.columns)
    }
}
