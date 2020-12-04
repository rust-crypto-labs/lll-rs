use rug;

use std::{cmp, ops};

pub(crate) trait Scalars {
    type Integer;
    type Fraction: cmp::PartialOrd<Self::Integer>
        + cmp::PartialOrd<Self::Fraction>
        + for<'a> ops::Mul<&'a Self::Fraction, Output = Self::Fraction>
        + for<'a> ops::Div<&'a Self::Fraction, Output = Self::Fraction>
        + for<'a> ops::SubAssign<&'a Self::Fraction>;

    fn from_frac(num: Self::Integer, denom: Self::Integer) -> Self::Fraction;
    fn from_int(n: Self::Integer) -> Self::Fraction;
    fn from_frac_i32(f: (i32, i32)) -> Self::Fraction;
    fn from_f64(f: f64) -> Option<Self::Fraction>;
    fn round(n: Self::Fraction) -> Self::Integer;
}

pub(crate) struct Float;
impl Scalars for Float {
    type Integer = f64;
    type Fraction = f64;

    fn from_frac(num: Self::Integer, denom: Self::Integer) -> Self::Fraction {
        num / denom
    }

    fn from_frac_i32(f: (i32, i32)) -> Self::Fraction {
        f64::from(f.0) / f64::from(f.1)
    }

    fn from_int(n: Self::Integer) -> Self::Fraction {
        n
    }

    fn round(f: Self::Fraction) -> Self::Integer {
        f.round()
    }

    fn from_f64(f: f64) -> Option<Self::Fraction> {
        if f.is_nan() || f.is_infinite() {
            None
        } else {
            Some(f)
        }
    }
}

pub(crate) struct BigNum;
impl Scalars for BigNum {
    type Integer = rug::Integer;
    type Fraction = rug::Rational;

    fn from_frac(num: Self::Integer, denom: Self::Integer) -> Self::Fraction {
        rug::Rational::from((num, denom))
    }

    fn from_frac_i32(f: (i32, i32)) -> Self::Fraction {
        rug::Rational::from(f)
    }

    fn from_int(n: Self::Integer) -> Self::Fraction {
        rug::Rational::from(n)
    }

    fn round(f: Self::Fraction) -> Self::Integer {
        f.fract_round(rug::Integer::new()).1
    }

    fn from_f64(f: f64) -> Option<Self::Fraction> {
        Self::Fraction::from_f64(f)
    }
}
