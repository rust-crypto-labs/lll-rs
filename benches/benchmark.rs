#[macro_use]
extern crate criterion;
extern crate lll_rs;
extern crate rug;

mod benchmarks {
    use criterion::Criterion;

    use lll_rs::{l2, lll, Matrix};

    pub fn bench_big_int_reduction_lll(c: &mut Criterion) {
        type I = rug::Integer;
        // "Bad" lattice basis
        let mut basis: Matrix<I> = Matrix::from_matrix(vec![
            vec![
                I::from(1) << 10000,
                I::from(0),
                I::from(0),
                I::from(1345) << 789,
            ],
            vec![
                I::from(0),
                I::from(1) << 500,
                I::from(0),
                I::from(35) << 3505,
            ],
            vec![
                I::from(0),
                I::from(0),
                I::from(1) << 1000,
                I::from(154) << 5000,
            ],
        ]);

        c.bench_function("lattice_reduce (biglll)", move |b| {
            b.iter(|| lll::lll_bignum(&mut basis))
        });
    }

    pub fn bench_big_int_reduction_l2(c: &mut Criterion) {
        type I = rug::Integer;
        // "Bad" lattice basis
        let mut basis: Matrix<I> = Matrix::from_matrix(vec![
            vec![
                I::from(1) << 10000,
                I::from(0),
                I::from(0),
                I::from(1345) << 789,
            ],
            vec![
                I::from(0),
                I::from(1) << 500,
                I::from(0),
                I::from(35) << 3505,
            ],
            vec![
                I::from(0),
                I::from(0),
                I::from(1) << 1000,
                I::from(154) << 5000,
            ],
        ]);

        c.bench_function("lattice_reduce (bigl2)", move |b| {
            b.iter(|| l2::lll_bignum(&mut basis, 0.501, 0.998))
        });
    }
}

criterion_group!(big_reduce_lll, benchmarks::bench_big_int_reduction_lll);
criterion_group!(big_reduce_l2, benchmarks::bench_big_int_reduction_l2);
criterion_main!(big_reduce_lll, big_reduce_l2);
