#[macro_use]
extern crate criterion;
extern crate lll_rs;
extern crate rug;

mod benchmarks {
    use criterion::Criterion;
    use lll_rs::lll;
    use lll_rs::l2;
    use lll_rs::matrix::Matrix;
    use lll_rs::vector::BigVector;
    use rug::Integer;

    pub fn bench_big_int_reduction_lll(c: &mut Criterion) {
        // "Bad" lattice basis
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0] = BigVector::from_vector(vec![
            Integer::from(1) << 10000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345) << 789,
        ]);
        basis[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1) << 500,
            Integer::from(0),
            Integer::from(35) << 3505,
        ]);
        basis[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1) << 1000,
            Integer::from(154) << 5000,
        ]);

        c.bench_function("big_lattice_reduce_lll", move |b| {
            b.iter(|| lll::big_lattice_reduce(&mut basis))
        });
    }

        pub fn bench_big_int_reduction_l2(c: &mut Criterion) {
        // "Bad" lattice basis
        let mut basis: Matrix<BigVector> = Matrix::init(3, 4);
        basis[0] = BigVector::from_vector(vec![
            Integer::from(1) << 10000,
            Integer::from(0),
            Integer::from(0),
            Integer::from(1345) << 789,
        ]);
        basis[1] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(1) << 500,
            Integer::from(0),
            Integer::from(35) << 3505,
        ]);
        basis[2] = BigVector::from_vector(vec![
            Integer::from(0),
            Integer::from(0),
            Integer::from(1) << 1000,
            Integer::from(154) << 5000,
        ]);

        c.bench_function("big_lattice_reduce_l2", move |b| {
            b.iter(|| l2::big_lattice_reduce(&mut basis, 0.501, 0.998))
        });
    }
}

criterion_group!(big_reduce_lll, benchmarks::bench_big_int_reduction_lll);
criterion_group!(big_reduce_l2, benchmarks::bench_big_int_reduction_l2);
criterion_main!(big_reduce_lll, big_reduce_l2);
