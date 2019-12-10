# lll-rs

`lll-rs` is an implementation of the Lenstra-Lenstra-Lovász lattice basis reduction algorithm (LLL [1]) in Rust as well as the implementation of an improved version, the L² algorithm (L² [2]). The library comes with a set of simple helpers to create vectors and matrices to perform lattice basis reduction.

## Theory

### Lattice

A lattice Λ is a dicrete subgroup of a vector space E, as an example if E = ℝⁿ (Lattice [3])

`X ∊ Λ <=> X = l_1 * b_1 + ... + l_n * b_n  with (l_i) in ℤ and (b_i) in ℝ`

[TODO]: Shortest vector

### LLL algorithm

Simply put, the LLL algorithm tries to solve the shortest vector problem by performing a Gram-Schimdt orthogonalization but at each step floors the value found to keep it in ℤ.
That way it can perform a base reduction algorithm in polynomial time

### Usage

LLL can be used for many operations:

- cracking RSA: The coppersmith attack (Coppersmith [4])
- Finding roots of polynomials with integer coefficients
- [TODO]: more

## The code

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

// Perfom the LLL basis redution
biglll::lattice_reduce(&mut basis);

// OR
// Perfom the LLL basis redution
// Specify the delta and eta coefficient for the reduction
bigl2::lattice_reduce(&mut basis, 0.5005, 0.999);
```

### Support

`lll-rs` supports:

- `rug::Integer` through `BigVector`
- `rug::Rational` through `RationalVector`
- `f64` through `VectorF` (simply performs a standard Gram-Schmidt orthogonalisation)

### Example

## References and documentation

- https://openaccess.leidenuniv.nl/bitstream/handle/1887/3810/346_050.pdf
- https://perso.ens-lyon.fr/damien.stehle/downloads/LLL25.pdf
- https://en.wikipedia.org/wiki/Lattice_(group)
- https://en.wikipedia.org/wiki/Coppersmith%27s_attack

[1]: https://openaccess.leidenuniv.nl/bitstream/handle/1887/3810/346_050.pdf
[2]: https://perso.ens-lyon.fr/damien.stehle/downloads/LLL25.pdf
[3]: https://en.wikipedia.org/wiki/Lattice_(group)
[4]: https://en.wikipedia.org/wiki/Coppersmith%27s_attack
