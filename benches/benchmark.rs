#[macro_use]
extern crate criterion;
extern crate lll_rs;
extern crate rug;

mod benchmarks {
    use criterion::Criterion;
    use lll_rs::lll::*;
    use lll_rs::matrix::*;
    use lll_rs::vector::*;
    use rug::*;

    pub fn bench_big_int_reduction(c: &mut Criterion) {
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

        c.bench_function("big_lattice_reduce", move |b| {
            b.iter(|| big_lattice_reduce(&mut basis))
        });
    }
}

criterion_group!(big_reduce, benchmarks::bench_big_int_reduction);
criterion_main!(big_reduce);
