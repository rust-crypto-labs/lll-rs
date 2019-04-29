//! A Rust implementation of the (basic) Lenstra-Lenstra-Lovasz lattice reduction algorithm
extern crate rug;

pub mod l2;
pub mod lll;
pub mod matrix;
pub mod vector;

#[cfg(test)]
mod test {
    use crate::l2;
    use crate::lll;
    use crate::matrix::Matrix;
    use crate::rug::Integer;
    use crate::vector::{BigVector, VectorF};

    #[test]
    fn test_lll() {
        // "Bad" lattice basis
        let mut basis: Matrix<VectorF> = Matrix::init(3, 4);
        basis[0] = VectorF::from_vector(vec![1., 0., 0., 1345.]);
        basis[1] = VectorF::from_vector(vec![0., 1., 0., 35.]);
        basis[2] = VectorF::from_vector(vec![0., 0., 1., 154.]);
        println!("{:?}", basis);

        // "Good" lattice basis
        lll::lattice_reduce(&mut basis);
        println!("{:?}", basis);
    }

    #[test]
    fn test_big_lll() {
        // "Bad" lattice basis
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0] = BigVector::from_vector(vec![
            Integer::from(1) << 100000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345),
        ]);
        basis[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(35),
        ]);
        basis[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(154),
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        lll::big_lattice_reduce(&mut basis);
        println!("{:?}", basis);
    }

    #[test]
    fn test_l2() {
        // "Bad" lattice basis
        let mut basis: Matrix<VectorF> = Matrix::init(3, 4);
        basis[0] = VectorF::from_vector(vec![1., 0., 0., 1345.]);
        basis[1] = VectorF::from_vector(vec![0., 1., 0., 35.]);
        basis[2] = VectorF::from_vector(vec![0., 0., 1., 154.]);
        println!("{:?}", basis);

        // "Good" lattice basis
        l2::lattice_reduce(&mut basis, 0.501, 0.998);
        println!("{:?}", basis);
    }

    #[test]
    fn test_big_l2() {
        // "Bad" lattice basis
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0] = BigVector::from_vector(vec![
            Integer::from(1) << 100000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345),
        ]);
        basis[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1),
            Integer::from(0),
            Integer::from(35),
        ]);
        basis[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1),
            Integer::from(154),
        ]);
        println!("{:?}", basis);

        // "Good" lattice basis
        l2::big_lattice_reduce(&mut basis, 0.5005, 0.999);
        println!("{:?}", basis);
    }
}
