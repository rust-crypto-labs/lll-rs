//! A Rust implementation of the (basic) Lenstra-Lenstra-Lovasz lattice reduction algorithm
extern crate rug;

pub mod l2;
pub mod lll;
pub mod matrix;
pub mod vector;

#[cfg(test)]
mod test {
    use crate::lll::{big_lattice_reduce, lattice_reduce};
    use crate::matrix::Matrix;
    use crate::rug::Integer;
    use crate::vector::{BigVector, VectorF};

    #[test]
    fn test_lll() {
        let matrix_type: Matrix<VectorF> = Matrix::init(3);

        // "Bad" lattice basis
        let mut basis = matrix_type.identity();
        basis.columns[0] = VectorF::from_vector(vec![1., 0., 0., 1345.]);
        basis.columns[1] = VectorF::from_vector(vec![0., 1., 0., 35.]);
        basis.columns[2] = VectorF::from_vector(vec![0., 0., 1., 154.]);
        println!("{:?}", basis);

        // "Good" lattice basis
        lattice_reduce(&mut basis);
        println!("{:?}", basis);
    }

    #[test]
    fn test_big_lll() {
        let matrix_type: Matrix<BigVector> = Matrix::init(3);

        // "Bad" lattice basis
        let mut basis = matrix_type.identity();
        basis.columns[0] = BigVector::from_vector(vec![
            Integer::from(1) << 100000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345),
        ]);
        basis.columns[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(35),
        ]);
        basis.columns[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(154),
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        big_lattice_reduce(&mut basis);
        println!("{:?}", basis);
    }
}
