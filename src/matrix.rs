use fmt::Debug;
use std::fmt;

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
}

impl<T> fmt::Debug for Matrix<T>
where
    T: Vector + Debug,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}\n", self.columns)
    }
}
