use rug::{Integer, Rational};
use std::{cmp, ops};

pub(crate) trait FromExt<T> {
    fn from_ext(_: T) -> Self;
}

macro_rules! impl_from_ext {
    ($from_type:ty, $to_type:ty, $code:expr) => {
        impl<'a> FromExt<$from_type> for $to_type {
            fn from_ext(f: $from_type) -> Self {
                $code(f)
            }
        }
    };
}

pub(crate) trait Scalars {
    type Integer;
    type Fraction: cmp::PartialOrd<Self::Integer>
        + for<'a> FromExt<&'a Self::Integer>
        + FromExt<f64>
        + FromExt<(Self::Integer, Self::Integer)>
        + FromExt<(i32, i32)>
        + cmp::PartialOrd<Self::Fraction>
        + for<'a> ops::Mul<&'a Self::Fraction, Output = Self::Fraction>
        + for<'a> ops::Div<&'a Self::Fraction, Output = Self::Fraction>
        + for<'a> ops::SubAssign<&'a Self::Fraction>;

    fn round(n: &Self::Fraction) -> Self::Integer;
    fn round_div(n: Self::Integer, d: Self::Integer) -> Self::Integer;
}

impl_from_ext!(&f64, f64, |f: &f64| *f);
impl_from_ext!((f64, f64), f64, |(n, d)| n / d);
impl_from_ext!(f64, f64, |f| f);
impl_from_ext!((i32, i32), f64, |(n, d)| f64::from(n) / f64::from(d));

pub(crate) struct Float;
impl Scalars for Float {
    type Integer = f64;
    type Fraction = f64;

    fn round(f: &Self::Fraction) -> Self::Integer {
        f.round()
    }

    fn round_div(n: Self::Integer, d: Self::Integer) -> Self::Integer {
        (n / d).round()
    }
}

impl_from_ext!(&Integer, Rational, |f: &Integer| Rational::from(f));
impl_from_ext!((Integer, Integer), Rational, |(n, d)| Rational::from((
    n, d
)));
impl_from_ext!(f64, Rational, |f: f64| Rational::from_f64(f).unwrap());
impl_from_ext!((i32, i32), Rational, |(n, d)| Rational::from((n, d)));

pub(crate) struct BigNum;
impl Scalars for BigNum {
    type Integer = rug::Integer;
    type Fraction = rug::Rational;

    fn round(f: &Self::Fraction) -> Self::Integer {
        f.round_ref().into()
    }

    fn round_div(mut n: Self::Integer, mut d: Self::Integer) -> Self::Integer {
        n.div_rem_round_mut(&mut d);
        n
    }
}
