# lll-rs

`lll-rs` is an implementation of the Lenstra–Lenstra–Lovász lattice basis reduction algorithm (LLL [[LLL82](#LLL82)], [1a], [1b]) in Rust.

## Supported algorithms

- LLL reduction [1a]
- L² reduction [2]
- Standard Gram-Schmidt orthogonalisation

The library comes with a set of simple helpers to create vectors and matrices, with the following entries:

- Integers (`BigVector`, relying on `rug::Integer`)
- Rationals (`RationalVector`, relying on `rug::Rational`)
- Small rationals (`VectorF`, relying on `f64`)

`lll-rs` is far from feature-complete and should be considered experimental. Users willing to use a stable and battle-tested library should
consider `fplll` instead [fplll].

## Lattice reduction

A lattice Λ is a dicrete subgroup of some vector space E. A typical example (see e.g. [3]) is E = ℝⁿ and

`X ∊ Λ <=> X = l_1 * b_1 + ... + l_n * b_n  with (l_i) in ℤ and (b_i) in ℝ`

Lattices are much studied mathematical structures on which we can formulate some useful problems [4]. Some of
these problems are simpler to solve when a "good basis" is known for the lattice. Conversely it is
difficult to solve them when only a "bad basis" is known.

Simply put, the LLL algorithm provides such a "good basis"; it roughly does so by performing a (variant of) rounded Gram-Schimdt orthogonalization on the "bad basis". 
Remarkably, this algorithm runs in polynomial time which makes it possible to solve several lattice problems efficiently.

Applications of LLL include:

- Cryptanalysis of lattice-based cryptosystems (e.g. NTRU)
- Cryptanalysis of pseudo-random number generators (e.g. LCG and truncated LCG)
- Cryptanalysis of RSA (e.g. Coppersmith's attack [5])
- Cryptanalysis of knapsack-based cryptosystems
- Finding mathematical counterexamples (e.g. Merten's conjecture)
- Finding roots of polynomials with integer coefficients
- Finding integer relations between constants
- Decoding of error correcting codes 

## Example

```rust
// Init the matrix with Integer
let mut basis: Matrix<BigVector> = Matrix::init(3, 4);

// Populate the matix
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

// Perfom the LLL basis reduction
biglll::lattice_reduce(&mut basis);

// OR
// Perfom the L2 basis reduction
// Specify the delta and eta coefficient for the reduction
bigl2::lattice_reduce(&mut basis, 0.5005, 0.999);
```

## References and documentation

<a name="LLL82">[LLL82]</a> A. K. Lenstra, H. W. Lenstra, Jr. and L. Lovasz. Factoring polynomials with rational coefficients. Math. Ann., 261: 515–534 (1982)

- https://openaccess.leidenuniv.nl/bitstream/handle/1887/3810/346_050.pdf
- https://en.wikipedia.org/wiki/Lenstra–Lenstra–Lovász_lattice_basis_reduction_algorithm
- https://perso.ens-lyon.fr/damien.stehle/downloads/LLL25.pdf
- https://en.wikipedia.org/wiki/Lattice_(group)
- https://en.wikipedia.org/wiki/Lattice_problem
- https://en.wikipedia.org/wiki/Coppersmith%27s_attack

[1a]: https://openaccess.leidenuniv.nl/bitstream/handle/1887/3810/346_050.pdf
[1b]: https://en.wikipedia.org/wikiLenstra–Lenstra–Lovász_lattice_basis_reduction_algorithm
[2]: https://perso.ens-lyon.fr/damien.stehle/downloads/LLL25.pdf
[3]: https://en.wikipedia.org/wiki/Lattice_(group)
[4]: https://en.wikipedia.org/wiki/Lattice_problem
[5]: https://en.wikipedia.org/wiki/Coppersmith%27s_attack
[fplll]: https://github.com/fplll/fplll