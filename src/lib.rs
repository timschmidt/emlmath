use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt;
use std::str::FromStr;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ComplexBall {
    pub mid: Complex,
    pub rad: f64,
}

impl Complex {
    pub const fn new(re: f64, im: f64) -> Self {
        Self { re, im }
    }

    pub const fn zero() -> Self {
        Self::new(0.0, 0.0)
    }

    pub const fn one() -> Self {
        Self::new(1.0, 0.0)
    }

    pub const fn i() -> Self {
        Self::new(0.0, 1.0)
    }

    pub fn is_zero(self) -> bool {
        self.re == 0.0 && self.im == 0.0
    }

    pub fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    pub fn arg(self) -> f64 {
        self.im.atan2(self.re)
    }

    pub fn conj(self) -> Self {
        Self::new(self.re, -self.im)
    }

    pub fn exp(self) -> Self {
        let scale = self.re.exp();
        if self.im == 0.0 {
            return Self::new(scale, 0.0);
        }
        Self::new(
            safe_scaled_component(scale, self.im.cos()),
            safe_scaled_component(scale, self.im.sin()),
        )
    }

    pub fn ln(self) -> Self {
        if self.is_zero() {
            return Self::new(f64::NEG_INFINITY, 0.0);
        }
        Self::new(self.abs().ln(), self.arg())
    }

    // EML needs a mirrored log branch so that the compiled ln(x) tree reproduces
    // the standard principal branch instead of flipping sign on the cut.
    pub fn eml_log(self) -> Self {
        if self.is_zero() {
            return Self::new(f64::NEG_INFINITY, 0.0);
        }
        Self::new(self.abs().ln(), -self.arg())
    }

    pub fn sqrt(self) -> Self {
        if self.is_zero() {
            return Self::zero();
        }
        let r = self.abs();
        let re = ((r + self.re) / 2.0).sqrt();
        let im_mag = ((r - self.re) / 2.0).sqrt();
        let im = if self.im < 0.0 { -im_mag } else { im_mag };
        Self::new(re, im)
    }

    pub fn sin(self) -> Self {
        Self::new(
            self.re.sin() * self.im.cosh(),
            self.re.cos() * self.im.sinh(),
        )
    }

    pub fn cos(self) -> Self {
        Self::new(
            self.re.cos() * self.im.cosh(),
            -self.re.sin() * self.im.sinh(),
        )
    }

    pub fn tan(self) -> Self {
        self.sin() / self.cos()
    }

    pub fn sinh(self) -> Self {
        Self::new(
            self.re.sinh() * self.im.cos(),
            self.re.cosh() * self.im.sin(),
        )
    }

    pub fn cosh(self) -> Self {
        Self::new(
            self.re.cosh() * self.im.cos(),
            self.re.sinh() * self.im.sin(),
        )
    }

    pub fn tanh(self) -> Self {
        self.sinh() / self.cosh()
    }

    pub fn powc(self, exponent: Self) -> Self {
        (exponent * self.ln()).exp()
    }

    pub fn asin(self) -> Self {
        let i = Self::i();
        -i * (i * self + (Self::one() - self * self).sqrt()).ln()
    }

    pub fn acos(self) -> Self {
        let i = Self::i();
        -i * (self + i * (Self::one() - self * self).sqrt()).ln()
    }

    pub fn atan(self) -> Self {
        let i = Self::i();
        let half_i = Self::new(0.0, 0.5);
        half_i * ((Self::one() - i * self).ln() - (Self::one() + i * self).ln())
    }

    pub fn asinh(self) -> Self {
        (self + (self * self + Self::one()).sqrt()).ln()
    }

    pub fn acosh(self) -> Self {
        (self + (self + Self::one()).sqrt() * (self - Self::one()).sqrt()).ln()
    }

    pub fn atanh(self) -> Self {
        (Self::one() / Self::new(2.0, 0.0))
            * ((Self::one() + self).ln() - (Self::one() - self).ln())
    }
}

impl ComplexBall {
    pub fn new(mid: Complex, rad: f64) -> Self {
        Self {
            mid,
            rad: rad.max(0.0),
        }
    }

    pub fn point(mid: Complex) -> Self {
        Self::new(mid, 0.0)
    }

    pub fn from_real(mid: f64, rad: f64) -> Self {
        Self::new(Complex::new(mid, 0.0), rad)
    }

    pub fn contains(self, value: Complex) -> bool {
        (value - self.mid).abs() <= self.rad
    }

    pub fn abs_upper(self) -> f64 {
        self.mid.abs() + self.rad
    }

    pub fn contains_zero(self) -> bool {
        self.mid.abs() <= self.rad
    }

    pub fn exp(self) -> Self {
        let mid = self.mid.exp();
        let rad = self.rad * (self.mid.re + self.rad).exp();
        Self::new(mid, rad)
    }

    pub fn ln(self) -> Result<Self, BallEvalError> {
        if self.contains_zero() {
            return Err(BallEvalError::ContainsSingularity("ln"));
        }
        if self.crosses_nonpositive_real_axis() {
            return Err(BallEvalError::CrossesBranchCut("ln"));
        }

        let clearance = self.mid.abs() - self.rad;
        let rad = self.rad / clearance;
        Ok(Self::new(self.mid.ln(), rad))
    }

    pub fn eml_log(self) -> Result<Self, BallEvalError> {
        if self.contains_zero() {
            return Err(BallEvalError::ContainsSingularity("eml_log"));
        }
        if self.crosses_nonpositive_real_axis() {
            return Err(BallEvalError::CrossesBranchCut("eml_log"));
        }

        let clearance = self.mid.abs() - self.rad;
        let rad = self.rad / clearance;
        Ok(Self::new(self.mid.eml_log(), rad))
    }

    pub fn reciprocal(self) -> Result<Self, BallEvalError> {
        if self.contains_zero() {
            return Err(BallEvalError::ContainsSingularity("reciprocal"));
        }

        let mid_abs = self.mid.abs();
        let clearance = mid_abs - self.rad;
        let rad = self.rad / (mid_abs * clearance);
        Ok(Self::new(Complex::one() / self.mid, rad))
    }

    pub fn checked_div(self, rhs: Self) -> Result<Self, BallEvalError> {
        Ok(self * rhs.reciprocal()?)
    }

    pub fn sqrt(self) -> Result<Self, BallEvalError> {
        let half = Self::point(Complex::new(0.5, 0.0));
        (half * self.ln()?).exp_checked()
    }

    pub fn sin(self) -> Result<Self, BallEvalError> {
        let j = -Self::i();
        let two_i = Self::point(Complex::new(2.0, 0.0)) * j;
        let ix = j * self;
        let neg_ix = -ix;
        ix.exp_checked()?.sub(neg_ix.exp_checked()?)?.checked_div(two_i)
    }

    pub fn cos(self) -> Result<Self, BallEvalError> {
        let j = -Self::i();
        let ix = j * self;
        let neg_ix = -ix;
        ix.exp_checked()?.add(neg_ix.exp_checked()?)?.checked_div(Self::point(Complex::new(2.0, 0.0)))
    }

    pub fn tan(self) -> Result<Self, BallEvalError> {
        self.sin()?.checked_div(self.cos()?)
    }

    pub fn sinh(self) -> Result<Self, BallEvalError> {
        self.exp_checked()?
            .sub((-self).exp_checked()?)?
            .checked_div(Self::point(Complex::new(2.0, 0.0)))
    }

    pub fn cosh(self) -> Result<Self, BallEvalError> {
        self.exp_checked()?
            .add((-self).exp_checked()?)?
            .checked_div(Self::point(Complex::new(2.0, 0.0)))
    }

    pub fn tanh(self) -> Result<Self, BallEvalError> {
        self.sinh()?.checked_div(self.cosh()?)
    }

    pub fn powc(self, exponent: Self) -> Result<Self, BallEvalError> {
        (exponent * self.ln()?).exp_checked()
    }

    pub fn asin(self) -> Result<Self, BallEvalError> {
        let i = Self::i();
        let one = Self::one();
        let inner = i * self + (one - self * self).sqrt()?;
        Ok((-i) * inner.ln()?)
    }

    pub fn acos(self) -> Result<Self, BallEvalError> {
        let i = Self::i();
        let one = Self::one();
        let inner = self + i * (one - self * self).sqrt()?;
        Ok((-i) * inner.ln()?)
    }

    pub fn atan(self) -> Result<Self, BallEvalError> {
        let i = Self::i();
        let half_i = Self::point(Complex::new(0.0, 0.5));
        let left = (Self::one() - i * self).ln()?;
        let right = (Self::one() + i * self).ln()?;
        Ok(half_i * (left - right))
    }

    pub fn asinh(self) -> Result<Self, BallEvalError> {
        (self + (self * self + Self::one()).sqrt()?).ln()
    }

    pub fn acosh(self) -> Result<Self, BallEvalError> {
        (self + (self + Self::one()).sqrt()? * (self - Self::one()).sqrt()?).ln()
    }

    pub fn atanh(self) -> Result<Self, BallEvalError> {
        let half = Self::point(Complex::new(0.5, 0.0));
        Ok(half * ((Self::one() + self).ln()? - (Self::one() - self).ln()?))
    }

    pub fn one() -> Self {
        Self::point(Complex::one())
    }

    pub fn i() -> Self {
        Self::point(Complex::i())
    }

    fn exp_checked(self) -> Result<Self, BallEvalError> {
        Ok(self.exp())
    }

    fn add(self, rhs: Self) -> Result<Self, BallEvalError> {
        Ok(self + rhs)
    }

    fn sub(self, rhs: Self) -> Result<Self, BallEvalError> {
        Ok(self - rhs)
    }

    fn crosses_nonpositive_real_axis(self) -> bool {
        distance_to_nonpositive_real_axis(self.mid) <= self.rad
    }
}

fn safe_scaled_component(scale: f64, factor: f64) -> f64 {
    if scale == 0.0 || factor == 0.0 {
        0.0
    } else {
        scale * factor
    }
}

impl std::ops::Add for Complex {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.re + rhs.re, self.im + rhs.im)
    }
}

impl std::ops::Sub for Complex {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.re - rhs.re, self.im - rhs.im)
    }
}

impl std::ops::Mul for Complex {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::new(
            self.re * rhs.re - self.im * rhs.im,
            self.re * rhs.im + self.im * rhs.re,
        )
    }
}

impl std::ops::Div for Complex {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        let denom = rhs.re * rhs.re + rhs.im * rhs.im;
        Self::new(
            (self.re * rhs.re + self.im * rhs.im) / denom,
            (self.im * rhs.re - self.re * rhs.im) / denom,
        )
    }
}

impl std::ops::Neg for Complex {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.re, -self.im)
    }
}

impl std::ops::Add for ComplexBall {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::new(self.mid + rhs.mid, self.rad + rhs.rad)
    }
}

impl std::ops::Sub for ComplexBall {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::new(self.mid - rhs.mid, self.rad + rhs.rad)
    }
}

impl std::ops::Mul for ComplexBall {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let rad = self.mid.abs() * rhs.rad + rhs.mid.abs() * self.rad + self.rad * rhs.rad;
        Self::new(self.mid * rhs.mid, rad)
    }
}

impl std::ops::Neg for ComplexBall {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new(-self.mid, self.rad)
    }
}

impl fmt::Display for Complex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.im == 0.0 {
            return write!(f, "{}", self.re);
        }
        if self.re == 0.0 {
            return write!(f, "{}i", self.im);
        }
        write!(f, "{} {:+}i", self.re, self.im)
    }
}

impl fmt::Display for ComplexBall {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} +/- {}", self.mid, self.rad)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EmlExpr {
    One,
    Var(String),
    Eml(Box<EmlExpr>, Box<EmlExpr>),
}

#[derive(Clone, Debug, PartialEq)]
pub enum ScientificExpr {
    One,
    Var(String),
    Exp(Box<ScientificExpr>),
    Ln(Box<ScientificExpr>),
    Neg(Box<ScientificExpr>),
    Add(Box<ScientificExpr>, Box<ScientificExpr>),
    Sub(Box<ScientificExpr>, Box<ScientificExpr>),
    Mul(Box<ScientificExpr>, Box<ScientificExpr>),
    Div(Box<ScientificExpr>, Box<ScientificExpr>),
    Pow(Box<ScientificExpr>, Box<ScientificExpr>),
    Sqrt(Box<ScientificExpr>),
    Sin(Box<ScientificExpr>),
    Cos(Box<ScientificExpr>),
    Tan(Box<ScientificExpr>),
    Sinh(Box<ScientificExpr>),
    Cosh(Box<ScientificExpr>),
    Tanh(Box<ScientificExpr>),
    Asin(Box<ScientificExpr>),
    Acos(Box<ScientificExpr>),
    Atan(Box<ScientificExpr>),
    Asinh(Box<ScientificExpr>),
    Acosh(Box<ScientificExpr>),
    Atanh(Box<ScientificExpr>),
}

impl ScientificExpr {
    pub fn one() -> Self {
        Self::One
    }

    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(name.into())
    }

    pub fn exp(expr: Self) -> Self {
        Self::Exp(Box::new(expr))
    }

    pub fn ln(expr: Self) -> Self {
        Self::Ln(Box::new(expr))
    }

    pub fn neg(expr: Self) -> Self {
        Self::Neg(Box::new(expr))
    }

    pub fn add(left: Self, right: Self) -> Self {
        Self::Add(Box::new(left), Box::new(right))
    }

    pub fn sub(left: Self, right: Self) -> Self {
        Self::Sub(Box::new(left), Box::new(right))
    }

    pub fn mul(left: Self, right: Self) -> Self {
        Self::Mul(Box::new(left), Box::new(right))
    }

    pub fn div(left: Self, right: Self) -> Self {
        Self::Div(Box::new(left), Box::new(right))
    }

    pub fn pow(base: Self, exponent: Self) -> Self {
        Self::Pow(Box::new(base), Box::new(exponent))
    }

    pub fn sqrt(expr: Self) -> Self {
        Self::Sqrt(Box::new(expr))
    }

    pub fn sin(expr: Self) -> Self {
        Self::Sin(Box::new(expr))
    }

    pub fn cos(expr: Self) -> Self {
        Self::Cos(Box::new(expr))
    }

    pub fn tan(expr: Self) -> Self {
        Self::Tan(Box::new(expr))
    }

    pub fn sinh(expr: Self) -> Self {
        Self::Sinh(Box::new(expr))
    }

    pub fn cosh(expr: Self) -> Self {
        Self::Cosh(Box::new(expr))
    }

    pub fn tanh(expr: Self) -> Self {
        Self::Tanh(Box::new(expr))
    }

    pub fn asin(expr: Self) -> Self {
        Self::Asin(Box::new(expr))
    }

    pub fn acos(expr: Self) -> Self {
        Self::Acos(Box::new(expr))
    }

    pub fn atan(expr: Self) -> Self {
        Self::Atan(Box::new(expr))
    }

    pub fn asinh(expr: Self) -> Self {
        Self::Asinh(Box::new(expr))
    }

    pub fn acosh(expr: Self) -> Self {
        Self::Acosh(Box::new(expr))
    }

    pub fn atanh(expr: Self) -> Self {
        Self::Atanh(Box::new(expr))
    }

    pub fn zero() -> Self {
        Self::ln(Self::one())
    }

    pub fn e() -> Self {
        Self::exp(Self::one())
    }

    pub fn two() -> Self {
        Self::add(Self::one(), Self::one())
    }

    pub fn half() -> Self {
        Self::div(Self::one(), Self::two())
    }

    pub fn minus_one() -> Self {
        Self::sub(Self::zero(), Self::one())
    }

    pub fn i() -> Self {
        Self::sqrt(Self::minus_one())
    }

    pub fn pi() -> Self {
        Self::mul(Self::neg(Self::i()), Self::ln(Self::minus_one()))
    }

    pub fn eval(&self, vars: &HashMap<String, Complex>) -> Result<Complex, EvalError> {
        Ok(match self {
            Self::One => Complex::one(),
            Self::Var(name) => vars
                .get(name)
                .copied()
                .ok_or_else(|| EvalError::UnknownVariable(name.clone()))?,
            Self::Exp(expr) => expr.eval(vars)?.exp(),
            Self::Ln(expr) => expr.eval(vars)?.ln(),
            Self::Neg(expr) => -expr.eval(vars)?,
            Self::Add(left, right) => left.eval(vars)? + right.eval(vars)?,
            Self::Sub(left, right) => left.eval(vars)? - right.eval(vars)?,
            Self::Mul(left, right) => left.eval(vars)? * right.eval(vars)?,
            Self::Div(left, right) => left.eval(vars)? / right.eval(vars)?,
            Self::Pow(base, exponent) => base.eval(vars)?.powc(exponent.eval(vars)?),
            Self::Sqrt(expr) => expr.eval(vars)?.sqrt(),
            Self::Sin(expr) => expr.eval(vars)?.sin(),
            Self::Cos(expr) => expr.eval(vars)?.cos(),
            Self::Tan(expr) => expr.eval(vars)?.tan(),
            Self::Sinh(expr) => expr.eval(vars)?.sinh(),
            Self::Cosh(expr) => expr.eval(vars)?.cosh(),
            Self::Tanh(expr) => expr.eval(vars)?.tanh(),
            Self::Asin(expr) => expr.eval(vars)?.asin(),
            Self::Acos(expr) => expr.eval(vars)?.acos(),
            Self::Atan(expr) => expr.eval(vars)?.atan(),
            Self::Asinh(expr) => expr.eval(vars)?.asinh(),
            Self::Acosh(expr) => expr.eval(vars)?.acosh(),
            Self::Atanh(expr) => expr.eval(vars)?.atanh(),
        })
    }

    pub fn eval_ball(
        &self,
        vars: &HashMap<String, ComplexBall>,
    ) -> Result<ComplexBall, BallEvalError> {
        Ok(match self {
            Self::One => ComplexBall::one(),
            Self::Var(name) => vars
                .get(name)
                .copied()
                .ok_or_else(|| BallEvalError::UnknownVariable(name.clone()))?,
            Self::Exp(expr) => expr.eval_ball(vars)?.exp(),
            Self::Ln(expr) => expr.eval_ball(vars)?.ln()?,
            Self::Neg(expr) => -expr.eval_ball(vars)?,
            Self::Add(left, right) => left.eval_ball(vars)? + right.eval_ball(vars)?,
            Self::Sub(left, right) => left.eval_ball(vars)? - right.eval_ball(vars)?,
            Self::Mul(left, right) => left.eval_ball(vars)? * right.eval_ball(vars)?,
            Self::Div(left, right) => left.eval_ball(vars)?.checked_div(right.eval_ball(vars)?)?,
            Self::Pow(base, exponent) => base.eval_ball(vars)?.powc(exponent.eval_ball(vars)?)?,
            Self::Sqrt(expr) => expr.eval_ball(vars)?.sqrt()?,
            Self::Sin(expr) => expr.eval_ball(vars)?.sin()?,
            Self::Cos(expr) => expr.eval_ball(vars)?.cos()?,
            Self::Tan(expr) => expr.eval_ball(vars)?.tan()?,
            Self::Sinh(expr) => expr.eval_ball(vars)?.sinh()?,
            Self::Cosh(expr) => expr.eval_ball(vars)?.cosh()?,
            Self::Tanh(expr) => expr.eval_ball(vars)?.tanh()?,
            Self::Asin(expr) => expr.eval_ball(vars)?.asin()?,
            Self::Acos(expr) => expr.eval_ball(vars)?.acos()?,
            Self::Atan(expr) => expr.eval_ball(vars)?.atan()?,
            Self::Asinh(expr) => expr.eval_ball(vars)?.asinh()?,
            Self::Acosh(expr) => expr.eval_ball(vars)?.acosh()?,
            Self::Atanh(expr) => expr.eval_ball(vars)?.atanh()?,
        })
    }

    pub fn to_eml(&self) -> EmlExpr {
        match self {
            Self::One => EmlExpr::one(),
            Self::Var(name) => EmlExpr::var(name.clone()),
            Self::Exp(expr) => EmlExpr::exp(expr.to_eml()),
            Self::Ln(expr) => EmlExpr::ln(expr.to_eml()),
            Self::Neg(expr) => EmlExpr::neg(expr.to_eml()),
            Self::Add(left, right) => EmlExpr::add(left.to_eml(), right.to_eml()),
            Self::Sub(left, right) => EmlExpr::sub(left.to_eml(), right.to_eml()),
            Self::Mul(left, right) => EmlExpr::mul(left.to_eml(), right.to_eml()),
            Self::Div(left, right) => EmlExpr::div(left.to_eml(), right.to_eml()),
            Self::Pow(base, exponent) => EmlExpr::pow(base.to_eml(), exponent.to_eml()),
            Self::Sqrt(expr) => EmlExpr::sqrt(expr.to_eml()),
            Self::Sin(expr) => EmlExpr::sin(expr.to_eml()),
            Self::Cos(expr) => EmlExpr::cos(expr.to_eml()),
            Self::Tan(expr) => EmlExpr::tan(expr.to_eml()),
            Self::Sinh(expr) => EmlExpr::sinh(expr.to_eml()),
            Self::Cosh(expr) => EmlExpr::cosh(expr.to_eml()),
            Self::Tanh(expr) => EmlExpr::tanh(expr.to_eml()),
            Self::Asin(expr) => EmlExpr::asin(expr.to_eml()),
            Self::Acos(expr) => EmlExpr::acos(expr.to_eml()),
            Self::Atan(expr) => EmlExpr::atan(expr.to_eml()),
            Self::Asinh(expr) => EmlExpr::asinh(expr.to_eml()),
            Self::Acosh(expr) => EmlExpr::acosh(expr.to_eml()),
            Self::Atanh(expr) => EmlExpr::atanh(expr.to_eml()),
        }
    }
}

impl EmlExpr {
    pub fn one() -> Self {
        Self::One
    }

    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(name.into())
    }

    pub fn eml(left: Self, right: Self) -> Self {
        Self::Eml(Box::new(left), Box::new(right))
    }

    pub fn eval(&self, vars: &HashMap<String, Complex>) -> Result<Complex, EvalError> {
        match self {
            Self::One => Ok(Complex::one()),
            Self::Var(name) => vars
                .get(name)
                .copied()
                .ok_or_else(|| EvalError::UnknownVariable(name.clone())),
            Self::Eml(left, right) => {
                let x = left.eval(vars)?;
                let y = right.eval(vars)?;
                Ok(x.exp() - y.eml_log())
            }
        }
    }

    pub fn eval_ball(
        &self,
        vars: &HashMap<String, ComplexBall>,
    ) -> Result<ComplexBall, BallEvalError> {
        match self {
            Self::One => Ok(ComplexBall::one()),
            Self::Var(name) => vars
                .get(name)
                .copied()
                .ok_or_else(|| BallEvalError::UnknownVariable(name.clone())),
            Self::Eml(left, right) => {
                let x = left.eval_ball(vars)?;
                let y = right.eval_ball(vars)?;
                Ok(x.exp() - y.eml_log()?)
            }
        }
    }

    pub fn node_count(&self) -> usize {
        match self {
            Self::One | Self::Var(_) => 1,
            Self::Eml(left, right) => 1 + left.node_count() + right.node_count(),
        }
    }

    pub fn exp(expr: Self) -> Self {
        Self::eml(expr, Self::one())
    }

    // Paper formula: ln(z) = eml(1, eml(eml(1, z), 1))
    pub fn ln(expr: Self) -> Self {
        Self::eml(
            Self::one(),
            Self::eml(Self::eml(Self::one(), expr), Self::one()),
        )
    }

    pub fn zero() -> Self {
        Self::ln(Self::one())
    }

    fn sub_raw(left: Self, right: Self) -> Self {
        Self::eml(Self::ln(left), Self::exp(right))
    }

    fn neg_raw(expr: Self) -> Self {
        Self::sub_raw(Self::zero(), expr)
    }

    pub fn sub(left: Self, right: Self) -> Self {
        Self::add(left, Self::neg(right))
    }

    pub fn neg(expr: Self) -> Self {
        Self::mul(Self::minus_one(), expr)
    }

    pub fn add(left: Self, right: Self) -> Self {
        Self::sub_raw(left, Self::neg_raw(right))
    }

    pub fn e() -> Self {
        Self::exp(Self::one())
    }

    pub fn reciprocal(expr: Self) -> Self {
        Self::exp(Self::neg(Self::ln(expr)))
    }

    pub fn mul(left: Self, right: Self) -> Self {
        Self::exp(Self::add(Self::ln(left), Self::ln(right)))
    }

    pub fn div(left: Self, right: Self) -> Self {
        Self::mul(left, Self::reciprocal(right))
    }

    pub fn pow(base: Self, exponent: Self) -> Self {
        Self::exp(Self::mul(exponent, Self::ln(base)))
    }

    pub fn two() -> Self {
        Self::add(Self::one(), Self::one())
    }

    pub fn half() -> Self {
        Self::reciprocal(Self::two())
    }

    pub fn minus_one() -> Self {
        Self::neg_raw(Self::one())
    }

    pub fn i() -> Self {
        Self::sqrt(Self::minus_one())
    }

    // Paper/compiler workaround: use the sign-corrected imaginary unit inside
    // derived EML formulas that rely on ln(-1), while keeping the public i()
    // constructor itself unchanged.
    fn compiler_i() -> Self {
        Self::neg(Self::i())
    }

    pub fn pi() -> Self {
        Self::mul(Self::neg(Self::i()), Self::ln(Self::minus_one()))
    }

    pub fn sqrt(expr: Self) -> Self {
        Self::pow(expr, Self::half())
    }

    pub fn sin(expr: Self) -> Self {
        let i = Self::compiler_i();
        let two_i = Self::mul(Self::two(), i.clone());
        let ix = Self::mul(i, expr.clone());
        let neg_ix = Self::neg(ix.clone());
        Self::div(Self::sub(Self::exp(ix), Self::exp(neg_ix)), two_i)
    }

    pub fn cos(expr: Self) -> Self {
        let i = Self::compiler_i();
        let ix = Self::mul(i, expr.clone());
        let neg_ix = Self::neg(ix.clone());
        Self::div(Self::add(Self::exp(ix), Self::exp(neg_ix)), Self::two())
    }

    pub fn tan(expr: Self) -> Self {
        Self::div(Self::sin(expr.clone()), Self::cos(expr))
    }

    pub fn sinh(expr: Self) -> Self {
        Self::div(
            Self::sub(Self::exp(expr.clone()), Self::exp(Self::neg(expr))),
            Self::two(),
        )
    }

    pub fn cosh(expr: Self) -> Self {
        Self::div(
            Self::add(Self::exp(expr.clone()), Self::exp(Self::neg(expr))),
            Self::two(),
        )
    }

    pub fn tanh(expr: Self) -> Self {
        Self::div(Self::sinh(expr.clone()), Self::cosh(expr))
    }

    pub fn asin(expr: Self) -> Self {
        let i = Self::compiler_i();
        let j = Self::neg(i.clone());
        let one = Self::one();
        let iz = Self::mul(j, expr.clone());
        let sqrt_term = Self::sqrt(Self::sub(one.clone(), Self::mul(expr.clone(), expr)));
        let inner = Self::add(iz, sqrt_term);
        Self::mul(i, Self::ln(inner))
    }

    pub fn acos(expr: Self) -> Self {
        let i = Self::compiler_i();
        let j = Self::neg(i.clone());
        let one = Self::one();
        let sqrt_term = Self::sqrt(Self::sub(one, Self::mul(expr.clone(), expr.clone())));
        let inner = Self::add(expr, Self::mul(j, sqrt_term));
        Self::mul(i, Self::ln(inner))
    }

    pub fn atan(expr: Self) -> Self {
        let i = Self::compiler_i();
        let j = Self::neg(i.clone());
        let one = Self::one();
        let half_i = Self::div(j.clone(), Self::two());
        let left = Self::ln(Self::sub(one.clone(), Self::mul(j.clone(), expr.clone())));
        let right = Self::ln(Self::add(one, Self::mul(j, expr)));
        Self::mul(half_i, Self::sub(left, right))
    }

    pub fn asinh(expr: Self) -> Self {
        let inner = Self::add(
            expr.clone(),
            Self::sqrt(Self::add(Self::mul(expr.clone(), expr), Self::one())),
        );
        Self::ln(inner)
    }

    pub fn acosh(expr: Self) -> Self {
        let inner = Self::add(
            expr.clone(),
            Self::mul(
                Self::sqrt(Self::add(expr.clone(), Self::one())),
                Self::sqrt(Self::sub(expr, Self::one())),
            ),
        );
        Self::ln(inner)
    }

    pub fn atanh(expr: Self) -> Self {
        Self::mul(
            Self::half(),
            Self::sub(
                Self::ln(Self::add(Self::one(), expr.clone())),
                Self::ln(Self::sub(Self::one(), expr)),
            ),
        )
    }
}

impl fmt::Display for EmlExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::One => write!(f, "1"),
            Self::Var(name) => write!(f, "{name}"),
            Self::Eml(left, right) => write!(f, "eml({left}, {right})"),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvalError {
    UnknownVariable(String),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BallEvalError {
    UnknownVariable(String),
    ContainsSingularity(&'static str),
    CrossesBranchCut(&'static str),
}

fn distance_to_nonpositive_real_axis(z: Complex) -> f64 {
    if z.re <= 0.0 { z.im.abs() } else { z.abs() }
}

pub fn eval_real(expr: &EmlExpr, assignments: &[(&str, f64)]) -> Result<Complex, EvalError> {
    let vars = assignments
        .iter()
        .map(|(name, value)| ((*name).to_string(), Complex::new(*value, 0.0)))
        .collect::<HashMap<_, _>>();
    expr.eval(&vars)
}

pub fn eval_real_scientific(
    expr: &ScientificExpr,
    assignments: &[(&str, f64)],
) -> Result<Complex, EvalError> {
    let vars = assignments
        .iter()
        .map(|(name, value)| ((*name).to_string(), Complex::new(*value, 0.0)))
        .collect::<HashMap<_, _>>();
    expr.eval(&vars)
}

pub fn eval_ball(
    expr: &EmlExpr,
    assignments: &[(&str, ComplexBall)],
) -> Result<ComplexBall, BallEvalError> {
    let vars = assignments
        .iter()
        .map(|(name, value)| ((*name).to_string(), *value))
        .collect::<HashMap<_, _>>();
    expr.eval_ball(&vars)
}

pub fn eval_ball_scientific(
    expr: &ScientificExpr,
    assignments: &[(&str, ComplexBall)],
) -> Result<ComplexBall, BallEvalError> {
    let vars = assignments
        .iter()
        .map(|(name, value)| ((*name).to_string(), *value))
        .collect::<HashMap<_, _>>();
    expr.eval_ball(&vars)
}

pub fn eval_real_ball(
    expr: &EmlExpr,
    assignments: &[(&str, f64, f64)],
) -> Result<ComplexBall, BallEvalError> {
    let vars = assignments
        .iter()
        .map(|(name, value, radius)| {
            ((*name).to_string(), ComplexBall::from_real(*value, *radius))
        })
        .collect::<HashMap<_, _>>();
    expr.eval_ball(&vars)
}

pub fn eval_real_ball_scientific(
    expr: &ScientificExpr,
    assignments: &[(&str, f64, f64)],
) -> Result<ComplexBall, BallEvalError> {
    let vars = assignments
        .iter()
        .map(|(name, value, radius)| {
            ((*name).to_string(), ComplexBall::from_real(*value, *radius))
        })
        .collect::<HashMap<_, _>>();
    expr.eval_ball(&vars)
}

pub fn x() -> EmlExpr {
    EmlExpr::var("x")
}

pub fn y() -> EmlExpr {
    EmlExpr::var("y")
}

pub fn z() -> EmlExpr {
    EmlExpr::var("z")
}

pub fn sx() -> ScientificExpr {
    ScientificExpr::var("x")
}

pub fn sy() -> ScientificExpr {
    ScientificExpr::var("y")
}

pub fn real_close(actual: Complex, expected: f64, tol: f64) -> bool {
    (actual.re - expected).abs() <= tol && actual.im.abs() <= tol
}

pub fn complex_close(actual: Complex, expected: Complex, tol: f64) -> bool {
    (actual.re - expected.re).abs() <= tol && (actual.im - expected.im).abs() <= tol
}

pub fn reference_pi() -> Complex {
    Complex::new(PI, 0.0)
}

pub fn reference_half_pi() -> Complex {
    Complex::new(FRAC_PI_2, 0.0)
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SurfaceScalarExpr {
    One,
    Var(String),
    E,
    I,
    Pi,
    Neg(Box<SurfaceScalarExpr>),
    Add(Box<SurfaceScalarExpr>, Box<SurfaceScalarExpr>),
    Sub(Box<SurfaceScalarExpr>, Box<SurfaceScalarExpr>),
    Mul(Box<SurfaceScalarExpr>, Box<SurfaceScalarExpr>),
    Div(Box<SurfaceScalarExpr>, Box<SurfaceScalarExpr>),
    Pow(Box<SurfaceScalarExpr>, Box<SurfaceScalarExpr>),
    Exp(Box<SurfaceScalarExpr>),
    Ln(Box<SurfaceScalarExpr>),
    Sqrt(Box<SurfaceScalarExpr>),
    Sin(Box<SurfaceScalarExpr>),
    Cos(Box<SurfaceScalarExpr>),
    Tan(Box<SurfaceScalarExpr>),
    Sinh(Box<SurfaceScalarExpr>),
    Cosh(Box<SurfaceScalarExpr>),
    Tanh(Box<SurfaceScalarExpr>),
    Asin(Box<SurfaceScalarExpr>),
    Acos(Box<SurfaceScalarExpr>),
    Atan(Box<SurfaceScalarExpr>),
    Asinh(Box<SurfaceScalarExpr>),
    Acosh(Box<SurfaceScalarExpr>),
    Atanh(Box<SurfaceScalarExpr>),
}

impl SurfaceScalarExpr {
    pub fn var(name: impl Into<String>) -> Self {
        Self::Var(name.into())
    }

    pub fn to_scientific(&self) -> ScientificExpr {
        match self {
            Self::One => ScientificExpr::one(),
            Self::Var(name) => ScientificExpr::var(name.clone()),
            Self::E => ScientificExpr::e(),
            Self::I => ScientificExpr::i(),
            Self::Pi => ScientificExpr::pi(),
            Self::Neg(expr) => ScientificExpr::neg(expr.to_scientific()),
            Self::Add(left, right) => {
                ScientificExpr::add(left.to_scientific(), right.to_scientific())
            }
            Self::Sub(left, right) => {
                ScientificExpr::sub(left.to_scientific(), right.to_scientific())
            }
            Self::Mul(left, right) => {
                ScientificExpr::mul(left.to_scientific(), right.to_scientific())
            }
            Self::Div(left, right) => {
                ScientificExpr::div(left.to_scientific(), right.to_scientific())
            }
            Self::Pow(base, exponent) => {
                ScientificExpr::pow(base.to_scientific(), exponent.to_scientific())
            }
            Self::Exp(expr) => ScientificExpr::exp(expr.to_scientific()),
            Self::Ln(expr) => ScientificExpr::ln(expr.to_scientific()),
            Self::Sqrt(expr) => ScientificExpr::sqrt(expr.to_scientific()),
            Self::Sin(expr) => ScientificExpr::sin(expr.to_scientific()),
            Self::Cos(expr) => ScientificExpr::cos(expr.to_scientific()),
            Self::Tan(expr) => ScientificExpr::tan(expr.to_scientific()),
            Self::Sinh(expr) => ScientificExpr::sinh(expr.to_scientific()),
            Self::Cosh(expr) => ScientificExpr::cosh(expr.to_scientific()),
            Self::Tanh(expr) => ScientificExpr::tanh(expr.to_scientific()),
            Self::Asin(expr) => ScientificExpr::asin(expr.to_scientific()),
            Self::Acos(expr) => ScientificExpr::acos(expr.to_scientific()),
            Self::Atan(expr) => ScientificExpr::atan(expr.to_scientific()),
            Self::Asinh(expr) => ScientificExpr::asinh(expr.to_scientific()),
            Self::Acosh(expr) => ScientificExpr::acosh(expr.to_scientific()),
            Self::Atanh(expr) => ScientificExpr::atanh(expr.to_scientific()),
        }
    }

    fn precedence(&self) -> u8 {
        match self {
            Self::Add(_, _) | Self::Sub(_, _) => 1,
            Self::Mul(_, _) | Self::Div(_, _) => 2,
            Self::Pow(_, _) => 3,
            Self::Neg(_) => 4,
            Self::Exp(_)
            | Self::Ln(_)
            | Self::Sqrt(_)
            | Self::Sin(_)
            | Self::Cos(_)
            | Self::Tan(_)
            | Self::Sinh(_)
            | Self::Cosh(_)
            | Self::Tanh(_)
            | Self::Asin(_)
            | Self::Acos(_)
            | Self::Atan(_)
            | Self::Asinh(_)
            | Self::Acosh(_)
            | Self::Atanh(_) => 5,
            Self::One | Self::Var(_) | Self::E | Self::I | Self::Pi => 6,
        }
    }

    fn fmt_with_prec(&self, f: &mut fmt::Formatter<'_>, parent_prec: u8) -> fmt::Result {
        let my_prec = self.precedence();
        let needs_parens = my_prec < parent_prec;
        if needs_parens {
            write!(f, "(")?;
        }
        match self {
            Self::One => write!(f, "1")?,
            Self::Var(name) => write!(f, "{name}")?,
            Self::E => write!(f, "e")?,
            Self::I => write!(f, "i")?,
            Self::Pi => write!(f, "pi")?,
            Self::Neg(expr) => {
                write!(f, "-")?;
                expr.fmt_with_prec(f, my_prec)?;
            }
            Self::Add(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " + ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::Sub(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " - ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::Mul(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " * ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::Div(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " / ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::Pow(base, exponent) => {
                base.fmt_with_prec(f, my_prec)?;
                write!(f, "^")?;
                exponent.fmt_with_prec(f, my_prec)?;
            }
            Self::Exp(expr) => write!(f, "exp({expr})")?,
            Self::Ln(expr) => write!(f, "ln({expr})")?,
            Self::Sqrt(expr) => write!(f, "sqrt({expr})")?,
            Self::Sin(expr) => write!(f, "sin({expr})")?,
            Self::Cos(expr) => write!(f, "cos({expr})")?,
            Self::Tan(expr) => write!(f, "tan({expr})")?,
            Self::Sinh(expr) => write!(f, "sinh({expr})")?,
            Self::Cosh(expr) => write!(f, "cosh({expr})")?,
            Self::Tanh(expr) => write!(f, "tanh({expr})")?,
            Self::Asin(expr) => write!(f, "asin({expr})")?,
            Self::Acos(expr) => write!(f, "acos({expr})")?,
            Self::Atan(expr) => write!(f, "atan({expr})")?,
            Self::Asinh(expr) => write!(f, "asinh({expr})")?,
            Self::Acosh(expr) => write!(f, "acosh({expr})")?,
            Self::Atanh(expr) => write!(f, "atanh({expr})")?,
        }
        if needs_parens {
            write!(f, ")")?;
        }
        Ok(())
    }
}

impl fmt::Display for SurfaceScalarExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.fmt_with_prec(f, 0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EmlCoordExpr {
    Scalar(SurfaceScalarExpr),
    Coord(SurfaceScalarExpr, SurfaceScalarExpr),
    Add(Box<EmlCoordExpr>, Box<EmlCoordExpr>),
    Sub(Box<EmlCoordExpr>, Box<EmlCoordExpr>),
    ScalarMul(Box<SurfaceScalarExpr>, Box<EmlCoordExpr>),
    Mul(Box<EmlCoordExpr>, Box<EmlCoordExpr>),
}

impl EmlCoordExpr {
    pub fn scalar(expr: SurfaceScalarExpr) -> Self {
        Self::Scalar(expr)
    }

    pub fn coord(a: SurfaceScalarExpr, b: SurfaceScalarExpr) -> Self {
        Self::Coord(a, b)
    }

    pub fn add(left: Self, right: Self) -> Self {
        Self::Add(Box::new(left), Box::new(right))
    }

    pub fn sub(left: Self, right: Self) -> Self {
        Self::Sub(Box::new(left), Box::new(right))
    }

    pub fn scalar_mul(scalar: SurfaceScalarExpr, expr: Self) -> Self {
        Self::ScalarMul(Box::new(scalar), Box::new(expr))
    }

    pub fn mul(left: Self, right: Self) -> Self {
        match (left.as_scalar(), right.as_scalar()) {
            (Some(left_scalar), Some(right_scalar)) => {
                Self::Scalar(SurfaceScalarExpr::Mul(
                    Box::new(left_scalar),
                    Box::new(right_scalar),
                ))
            }
            (Some(left_scalar), None) => Self::scalar_mul(left_scalar, right),
            (None, Some(right_scalar)) => Self::scalar_mul(right_scalar, left),
            (None, None) => Self::Mul(Box::new(left), Box::new(right)),
        }
    }

    pub fn to_scientific(&self) -> ScientificExpr {
        match self {
            Self::Scalar(expr) => expr.to_scientific(),
            Self::Coord(a, b) => {
                ScientificExpr::sub(ScientificExpr::exp(a.to_scientific()), ScientificExpr::ln(b.to_scientific()))
            }
            Self::Add(left, right) => {
                ScientificExpr::add(left.to_scientific(), right.to_scientific())
            }
            Self::Sub(left, right) => {
                ScientificExpr::sub(left.to_scientific(), right.to_scientific())
            }
            Self::ScalarMul(scalar, expr) => {
                ScientificExpr::mul(scalar.to_scientific(), expr.to_scientific())
            }
            Self::Mul(left, right) => {
                ScientificExpr::mul(left.to_scientific(), right.to_scientific())
            }
        }
    }

    pub fn to_eml(&self) -> EmlExpr {
        self.to_scientific().to_eml()
    }

    fn as_scalar(&self) -> Option<SurfaceScalarExpr> {
        match self {
            Self::Scalar(expr) => Some(expr.clone()),
            _ => None,
        }
    }

    fn precedence(&self) -> u8 {
        match self {
            Self::Add(_, _) | Self::Sub(_, _) => 1,
            Self::ScalarMul(_, _) | Self::Mul(_, _) => 2,
            Self::Coord(_, _) => 1,
            Self::Scalar(_) => 3,
        }
    }

    fn fmt_with_prec(&self, f: &mut fmt::Formatter<'_>, parent_prec: u8) -> fmt::Result {
        let my_prec = self.precedence();
        let needs_parens = my_prec < parent_prec;
        if needs_parens {
            write!(f, "(")?;
        }
        match self {
            Self::Scalar(expr) => expr.fmt_with_prec(f, 0)?,
            Self::Coord(a, b) => {
                a.fmt_with_prec(f, 0)?;
                match b {
                    SurfaceScalarExpr::Neg(inner) => {
                        write!(f, " - ε ")?;
                        inner.fmt_with_prec(f, 0)?;
                    }
                    _ => {
                        write!(f, " + ε ")?;
                        b.fmt_with_prec(f, 0)?;
                    }
                }
            }
            Self::Add(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " + ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::Sub(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " - ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::ScalarMul(scalar, expr) => {
                scalar.fmt_with_prec(f, my_prec)?;
                write!(f, " * ")?;
                expr.fmt_with_prec(f, my_prec + 1)?;
            }
            Self::Mul(left, right) => {
                left.fmt_with_prec(f, my_prec)?;
                write!(f, " * ")?;
                right.fmt_with_prec(f, my_prec + 1)?;
            }
        }
        if needs_parens {
            write!(f, ")")?;
        }
        Ok(())
    }
}

impl fmt::Display for EmlCoordExpr {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.fmt_with_prec(f, 0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SurfaceParseError {
    pub message: String,
}

impl SurfaceParseError {
    fn new(message: impl Into<String>) -> Self {
        Self {
            message: message.into(),
        }
    }
}

impl fmt::Display for SurfaceParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for SurfaceParseError {}

#[derive(Clone, Debug, PartialEq, Eq)]
enum SurfaceToken {
    One,
    Ident(String),
    Epsilon,
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    LParen,
    RParen,
}

fn tokenize_surface(input: &str) -> Result<Vec<SurfaceToken>, SurfaceParseError> {
    let mut tokens = Vec::new();
    let mut chars = input.chars().peekable();

    while let Some(ch) = chars.peek().copied() {
        match ch {
            ' ' | '\t' | '\n' | '\r' => {
                chars.next();
            }
            '1' => {
                chars.next();
                tokens.push(SurfaceToken::One);
            }
            '+' => {
                chars.next();
                tokens.push(SurfaceToken::Plus);
            }
            '-' => {
                chars.next();
                tokens.push(SurfaceToken::Minus);
            }
            '*' => {
                chars.next();
                tokens.push(SurfaceToken::Star);
            }
            '/' => {
                chars.next();
                tokens.push(SurfaceToken::Slash);
            }
            '^' => {
                chars.next();
                tokens.push(SurfaceToken::Caret);
            }
            '(' => {
                chars.next();
                tokens.push(SurfaceToken::LParen);
            }
            ')' => {
                chars.next();
                tokens.push(SurfaceToken::RParen);
            }
            'ε' => {
                chars.next();
                tokens.push(SurfaceToken::Epsilon);
            }
            c if c.is_ascii_alphabetic() || c == '_' => {
                let mut ident = String::new();
                while let Some(next) = chars.peek().copied() {
                    if next.is_ascii_alphanumeric() || next == '_' {
                        ident.push(next);
                        chars.next();
                    } else {
                        break;
                    }
                }
                match ident.as_str() {
                    "eps" | "epsilon" => tokens.push(SurfaceToken::Epsilon),
                    _ => tokens.push(SurfaceToken::Ident(ident)),
                }
            }
            _ => {
                return Err(SurfaceParseError::new(format!(
                    "unexpected character `{ch}` in surface expression"
                )));
            }
        }
    }

    Ok(tokens)
}

struct SurfaceParser {
    tokens: Vec<SurfaceToken>,
    pos: usize,
}

impl SurfaceParser {
    fn new(tokens: Vec<SurfaceToken>) -> Self {
        Self { tokens, pos: 0 }
    }

    fn parse_coord_expr(&mut self) -> Result<EmlCoordExpr, SurfaceParseError> {
        let expr = self.parse_coord_sum()?;
        if self.pos != self.tokens.len() {
            return Err(SurfaceParseError::new("unexpected trailing tokens"));
        }
        Ok(expr)
    }

    fn parse_coord_sum(&mut self) -> Result<EmlCoordExpr, SurfaceParseError> {
        let mut expr = self.parse_coord_product()?;
        loop {
            match self.peek() {
                Some(SurfaceToken::Plus) => {
                    self.pos += 1;
                    let rhs = self.parse_coord_product()?;
                    expr = EmlCoordExpr::add(expr, rhs);
                }
                Some(SurfaceToken::Minus) => {
                    self.pos += 1;
                    let rhs = self.parse_coord_product()?;
                    expr = EmlCoordExpr::sub(expr, rhs);
                }
                _ => return Ok(expr),
            }
        }
    }

    fn parse_coord_product(&mut self) -> Result<EmlCoordExpr, SurfaceParseError> {
        let mut expr = self.parse_coord_atom()?;
        while matches!(self.peek(), Some(SurfaceToken::Star)) {
            self.pos += 1;
            let rhs = self.parse_coord_atom()?;
            expr = EmlCoordExpr::mul(expr, rhs);
        }
        Ok(expr)
    }

    fn parse_coord_atom(&mut self) -> Result<EmlCoordExpr, SurfaceParseError> {
        if matches!(self.peek(), Some(SurfaceToken::LParen)) {
            self.pos += 1;
            let inner = self.parse_coord_sum()?;
            self.expect(SurfaceToken::RParen)?;
            return Ok(inner);
        }

        let scalar = self.parse_scalar_expr(true, true)?;
        let coord = match (self.peek(), self.peek_n(1)) {
            (Some(SurfaceToken::Plus), Some(SurfaceToken::Epsilon)) => {
                self.pos += 2;
                let rhs = self.parse_scalar_expr(false, false)?;
                EmlCoordExpr::coord(scalar, rhs)
            }
            (Some(SurfaceToken::Minus), Some(SurfaceToken::Epsilon)) => {
                self.pos += 2;
                let rhs = self.parse_scalar_expr(false, false)?;
                EmlCoordExpr::coord(scalar, SurfaceScalarExpr::Neg(Box::new(rhs)))
            }
            _ => EmlCoordExpr::scalar(scalar),
        };
        Ok(coord)
    }

    fn parse_scalar_expr(
        &mut self,
        stop_before_coord_suffix: bool,
        stop_before_coord_product: bool,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        self.parse_scalar_add_sub(stop_before_coord_suffix, stop_before_coord_product)
    }

    fn parse_scalar_add_sub(
        &mut self,
        stop_before_coord_suffix: bool,
        stop_before_coord_product: bool,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        let mut expr =
            self.parse_scalar_mul_div(stop_before_coord_suffix, stop_before_coord_product)?;
        loop {
            match self.peek() {
                Some(SurfaceToken::Plus)
                    if stop_before_coord_suffix
                        && matches!(self.peek_n(1), Some(SurfaceToken::Epsilon)) =>
                {
                    return Ok(expr);
                }
                Some(SurfaceToken::Minus)
                    if stop_before_coord_suffix
                        && matches!(self.peek_n(1), Some(SurfaceToken::Epsilon)) =>
                {
                    return Ok(expr);
                }
                Some(SurfaceToken::Plus) => {
                    self.pos += 1;
                    let rhs = self
                        .parse_scalar_mul_div(stop_before_coord_suffix, stop_before_coord_product)?;
                    expr = SurfaceScalarExpr::Add(Box::new(expr), Box::new(rhs));
                }
                Some(SurfaceToken::Minus) => {
                    self.pos += 1;
                    let rhs = self
                        .parse_scalar_mul_div(stop_before_coord_suffix, stop_before_coord_product)?;
                    expr = SurfaceScalarExpr::Sub(Box::new(expr), Box::new(rhs));
                }
                _ => return Ok(expr),
            }
        }
    }

    fn parse_scalar_mul_div(
        &mut self,
        stop_before_coord_suffix: bool,
        stop_before_coord_product: bool,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        let mut expr = self.parse_scalar_pow(stop_before_coord_suffix, stop_before_coord_product)?;
        loop {
            match self.peek() {
                Some(SurfaceToken::Star) if stop_before_coord_product => return Ok(expr),
                Some(SurfaceToken::Star) => {
                    self.pos += 1;
                    let rhs =
                        self.parse_scalar_pow(stop_before_coord_suffix, stop_before_coord_product)?;
                    expr = SurfaceScalarExpr::Mul(Box::new(expr), Box::new(rhs));
                }
                Some(SurfaceToken::Slash) => {
                    self.pos += 1;
                    let rhs =
                        self.parse_scalar_pow(stop_before_coord_suffix, stop_before_coord_product)?;
                    expr = SurfaceScalarExpr::Div(Box::new(expr), Box::new(rhs));
                }
                _ => return Ok(expr),
            }
        }
    }

    fn parse_scalar_pow(
        &mut self,
        stop_before_coord_suffix: bool,
        stop_before_coord_product: bool,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        let left =
            self.parse_scalar_unary(stop_before_coord_suffix, stop_before_coord_product)?;
        if matches!(self.peek(), Some(SurfaceToken::Caret)) {
            self.pos += 1;
            let right =
                self.parse_scalar_pow(stop_before_coord_suffix, stop_before_coord_product)?;
            return Ok(SurfaceScalarExpr::Pow(Box::new(left), Box::new(right)));
        }
        Ok(left)
    }

    fn parse_scalar_unary(
        &mut self,
        stop_before_coord_suffix: bool,
        stop_before_coord_product: bool,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        if matches!(self.peek(), Some(SurfaceToken::Minus)) {
            self.pos += 1;
            return Ok(SurfaceScalarExpr::Neg(Box::new(
                self.parse_scalar_unary(stop_before_coord_suffix, stop_before_coord_product)?,
            )));
        }
        self.parse_scalar_primary(stop_before_coord_suffix, stop_before_coord_product)
    }

    fn parse_scalar_primary(
        &mut self,
        stop_before_coord_suffix: bool,
        stop_before_coord_product: bool,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        match self.next() {
            Some(SurfaceToken::One) => Ok(SurfaceScalarExpr::One),
            Some(SurfaceToken::Ident(name)) => {
                if matches!(self.peek(), Some(SurfaceToken::LParen)) {
                    self.pos += 1;
                    let arg = self.parse_scalar_expr(false, false)?;
                    self.expect(SurfaceToken::RParen)?;
                    self.build_function(name, arg)
                } else {
                    match name.as_str() {
                        "e" => Ok(SurfaceScalarExpr::E),
                        "i" => Ok(SurfaceScalarExpr::I),
                        "pi" => Ok(SurfaceScalarExpr::Pi),
                        _ => Ok(SurfaceScalarExpr::Var(name)),
                    }
                }
            }
            Some(SurfaceToken::LParen) => {
                let expr = self.parse_scalar_expr(
                    stop_before_coord_suffix,
                    stop_before_coord_product,
                )?;
                self.expect(SurfaceToken::RParen)?;
                Ok(expr)
            }
            Some(SurfaceToken::Epsilon) => {
                Err(SurfaceParseError::new("ε must appear as `a + ε b` or `a - ε b`"))
            }
            other => Err(SurfaceParseError::new(format!(
                "unexpected token {:?} while parsing scalar expression",
                other
            ))),
        }
    }

    fn build_function(
        &self,
        name: String,
        arg: SurfaceScalarExpr,
    ) -> Result<SurfaceScalarExpr, SurfaceParseError> {
        let boxed = Box::new(arg);
        match name.as_str() {
            "exp" => Ok(SurfaceScalarExpr::Exp(boxed)),
            "ln" => Ok(SurfaceScalarExpr::Ln(boxed)),
            "sqrt" => Ok(SurfaceScalarExpr::Sqrt(boxed)),
            "sin" => Ok(SurfaceScalarExpr::Sin(boxed)),
            "cos" => Ok(SurfaceScalarExpr::Cos(boxed)),
            "tan" => Ok(SurfaceScalarExpr::Tan(boxed)),
            "sinh" => Ok(SurfaceScalarExpr::Sinh(boxed)),
            "cosh" => Ok(SurfaceScalarExpr::Cosh(boxed)),
            "tanh" => Ok(SurfaceScalarExpr::Tanh(boxed)),
            "asin" => Ok(SurfaceScalarExpr::Asin(boxed)),
            "acos" => Ok(SurfaceScalarExpr::Acos(boxed)),
            "atan" => Ok(SurfaceScalarExpr::Atan(boxed)),
            "asinh" => Ok(SurfaceScalarExpr::Asinh(boxed)),
            "acosh" => Ok(SurfaceScalarExpr::Acosh(boxed)),
            "atanh" => Ok(SurfaceScalarExpr::Atanh(boxed)),
            _ => Err(SurfaceParseError::new(format!(
                "unknown scalar function `{name}`"
            ))),
        }
    }

    fn expect(&mut self, expected: SurfaceToken) -> Result<(), SurfaceParseError> {
        let next = self.next();
        if next == Some(expected.clone()) {
            Ok(())
        } else {
            Err(SurfaceParseError::new(format!(
                "expected {:?}, found {:?}",
                expected, next
            )))
        }
    }

    fn next(&mut self) -> Option<SurfaceToken> {
        let token = self.tokens.get(self.pos).cloned();
        if token.is_some() {
            self.pos += 1;
        }
        token
    }

    fn peek(&self) -> Option<&SurfaceToken> {
        self.tokens.get(self.pos)
    }

    fn peek_n(&self, offset: usize) -> Option<&SurfaceToken> {
        self.tokens.get(self.pos + offset)
    }
}

impl FromStr for EmlCoordExpr {
    type Err = SurfaceParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tokens = tokenize_surface(s)?;
        SurfaceParser::new(tokens).parse_coord_expr()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_close(actual: Complex, expected: Complex, tol: f64) {
        assert!(
            complex_close(actual, expected, tol),
            "actual={actual:?} expected={expected:?} tol={tol}"
        );
    }

    fn assert_ball_contains(ball: ComplexBall, expected: Complex) {
        assert!(
            ball.contains(expected),
            "ball={ball:?} expected point={expected:?}"
        );
    }

    fn assert_eml_matches_scientific(
        scientific: ScientificExpr,
        vars: &[(&str, Complex)],
        tol: f64,
    ) {
        let vars = vars
            .iter()
            .map(|(name, value)| ((*name).to_string(), *value))
            .collect::<HashMap<_, _>>();
        let scientific_value = scientific.eval(&vars).unwrap();
        let eml_value = scientific.to_eml().eval(&vars).unwrap();
        assert!(
            complex_close(eml_value, scientific_value, tol),
            "scientific={scientific:?} eml={} scientific_value={scientific_value:?} eml_value={eml_value:?} tol={tol}",
            scientific.to_eml(),
        );
    }

    #[test]
    fn constants_work() {
        assert_close(
            eval_real(&EmlExpr::zero(), &[]).unwrap(),
            Complex::zero(),
            1e-12,
        );
        assert_close(
            eval_real(&EmlExpr::e(), &[]).unwrap(),
            Complex::new(std::f64::consts::E, 0.0),
            1e-12,
        );
        assert_close(
            eval_real(&EmlExpr::two(), &[]).unwrap(),
            Complex::new(2.0, 0.0),
            1e-10,
        );
        assert_close(
            eval_real(&EmlExpr::minus_one(), &[]).unwrap(),
            Complex::new(-1.0, 0.0),
            1e-9,
        );
        assert_close(
            eval_real(&EmlExpr::half(), &[]).unwrap(),
            Complex::new(0.5, 0.0),
            1e-9,
        );
    }

    #[test]
    fn arithmetic_works() {
        assert_close(
            eval_real(&EmlExpr::add(x(), y()), &[("x", 2.5), ("y", -0.75)]).unwrap(),
            Complex::new(1.75, 0.0),
            1e-9,
        );
        assert_close(
            eval_real(&EmlExpr::sub(x(), y()), &[("x", 2.5), ("y", -0.75)]).unwrap(),
            Complex::new(3.25, 0.0),
            1e-9,
        );
        assert_close(
            eval_real(&EmlExpr::mul(x(), y()), &[("x", 2.5), ("y", -0.75)]).unwrap(),
            Complex::new(-1.875, 0.0),
            1e-8,
        );
        assert_close(
            eval_real(&EmlExpr::div(x(), y()), &[("x", 2.5), ("y", -0.75)]).unwrap(),
            Complex::new(-3.3333333333333335, 0.0),
            1e-8,
        );
        assert_close(
            eval_real(&EmlExpr::pow(x(), y()), &[("x", 2.5), ("y", -0.75)]).unwrap(),
            Complex::new(2.5_f64.powf(-0.75), 0.0),
            1e-8,
        );
        assert_close(
            eval_real(&EmlExpr::sqrt(x()), &[("x", 2.5)]).unwrap(),
            Complex::new(2.5_f64.sqrt(), 0.0),
            1e-8,
        );
    }

    #[test]
    fn trig_works() {
        let x_value = Complex::new(0.4, 0.0);
        let vars = HashMap::from([("x".to_string(), x_value)]);

        assert_close(
            ScientificExpr::sin(sx()).eval(&vars).unwrap(),
            x_value.sin(),
            1e-12,
        );
        assert_close(
            ScientificExpr::cos(sx()).eval(&vars).unwrap(),
            x_value.cos(),
            1e-12,
        );
        assert_close(
            ScientificExpr::tan(sx()).eval(&vars).unwrap(),
            x_value.tan(),
            1e-12,
        );
        assert_close(
            ScientificExpr::sinh(sx()).eval(&vars).unwrap(),
            x_value.sinh(),
            1e-12,
        );
        assert_close(
            ScientificExpr::cosh(sx()).eval(&vars).unwrap(),
            x_value.cosh(),
            1e-12,
        );
        assert_close(
            ScientificExpr::tanh(sx()).eval(&vars).unwrap(),
            x_value.tanh(),
            1e-12,
        );
        assert!(ScientificExpr::sin(sx()).to_eml().node_count() > 1);
        assert!(ScientificExpr::tan(sx()).to_eml().node_count() > 1);
    }

    #[test]
    fn inverse_functions_work() {
        let x_value = Complex::new(0.25, 0.0);
        let vars = HashMap::from([("x".to_string(), x_value)]);

        assert_close(
            ScientificExpr::ln(sx()).eval(&vars).unwrap(),
            x_value.ln(),
            1e-12,
        );
        assert_close(
            ScientificExpr::asin(sx()).eval(&vars).unwrap(),
            x_value.asin(),
            1e-12,
        );
        assert_close(
            ScientificExpr::acos(sx()).eval(&vars).unwrap(),
            x_value.acos(),
            1e-12,
        );
        assert_close(
            ScientificExpr::atan(sx()).eval(&vars).unwrap(),
            x_value.atan(),
            1e-12,
        );
        assert_close(
            ScientificExpr::asinh(sx()).eval(&vars).unwrap(),
            x_value.asinh(),
            1e-12,
        );
        assert_close(
            ScientificExpr::acosh(ScientificExpr::add(sx(), ScientificExpr::two()))
                .eval(&vars)
                .unwrap(),
            (x_value + Complex::new(2.0, 0.0)).acosh(),
            1e-12,
        );
        assert_close(
            ScientificExpr::atanh(sx()).eval(&vars).unwrap(),
            x_value.atanh(),
            1e-12,
        );
        assert!(ScientificExpr::asin(sx()).to_eml().node_count() > 1);
    }

    #[test]
    fn derived_pi_and_i_work() {
        assert_close(
            eval_real_scientific(&ScientificExpr::i(), &[]).unwrap(),
            Complex::i(),
            1e-12,
        );
        assert_close(
            eval_real_scientific(&ScientificExpr::pi(), &[]).unwrap(),
            Complex::new(PI, 0.0),
            1e-12,
        );
        assert!(ScientificExpr::pi().to_eml().node_count() > 1);
    }

    #[test]
    fn compiled_constants_match_scientific() {
        for scientific in [
            ScientificExpr::zero(),
            ScientificExpr::e(),
            ScientificExpr::two(),
            ScientificExpr::half(),
            ScientificExpr::minus_one(),
            ScientificExpr::i(),
            ScientificExpr::pi(),
        ] {
            assert_eml_matches_scientific(scientific, &[], 1e-9);
        }
    }

    #[test]
    fn compiled_arithmetic_matches_scientific() {
        let vars = [
            ("x", Complex::new(2.5, 0.0)),
            ("y", Complex::new(-0.75, 0.0)),
        ];

        for scientific in [
            ScientificExpr::add(sx(), sy()),
            ScientificExpr::sub(sx(), sy()),
            ScientificExpr::mul(sx(), sy()),
            ScientificExpr::div(sx(), sy()),
            ScientificExpr::pow(sx(), sy()),
            ScientificExpr::sqrt(sx()),
        ] {
            assert_eml_matches_scientific(scientific, &vars, 1e-8);
        }
    }

    #[test]
    fn compiled_trig_and_hyperbolic_match_scientific() {
        let vars = [("x", Complex::new(0.4, 0.0))];

        for scientific in [
            ScientificExpr::sin(sx()),
            ScientificExpr::cos(sx()),
            ScientificExpr::tan(sx()),
            ScientificExpr::sinh(sx()),
            ScientificExpr::cosh(sx()),
            ScientificExpr::tanh(sx()),
        ] {
            assert_eml_matches_scientific(scientific, &vars, 1e-8);
        }
    }

    #[test]
    fn compiled_inverse_functions_match_scientific() {
        let vars = [("x", Complex::new(0.25, 0.0))];

        for scientific in [
            ScientificExpr::ln(sx()),
            ScientificExpr::asin(sx()),
            ScientificExpr::acos(sx()),
            ScientificExpr::atan(sx()),
            ScientificExpr::asinh(sx()),
            ScientificExpr::acosh(ScientificExpr::add(sx(), ScientificExpr::two())),
            ScientificExpr::atanh(sx()),
        ] {
            assert_eml_matches_scientific(scientific, &vars, 1e-8);
        }
    }

    #[test]
    fn mirrored_log_branch_matches_expected_signs() {
        let vars_upper = HashMap::from([("z".to_string(), Complex::new(-1.0, 0.5))]);
        let vars_lower = HashMap::from([("z".to_string(), Complex::new(-1.0, -0.5))]);
        let expr = EmlExpr::eml(EmlExpr::zero(), EmlExpr::var("z"));

        let upper = expr.eval(&vars_upper).unwrap();
        let lower = expr.eval(&vars_lower).unwrap();

        assert_close(
            upper,
            Complex::one() - Complex::new(-1.0, 0.5).eml_log(),
            1e-12,
        );
        assert_close(
            lower,
            Complex::one() - Complex::new(-1.0, -0.5).eml_log(),
            1e-12,
        );
        assert!(upper.im > 0.0);
        assert!(lower.im < 0.0);
    }

    #[test]
    fn compiled_ln_matches_scientific_across_branch_sides() {
        for value in [Complex::new(-1.0, 0.5), Complex::new(-1.0, -0.5)] {
            assert_eml_matches_scientific(ScientificExpr::ln(sx()), &[("x", value)], 1e-8);
        }
    }

    #[test]
    fn coord_surface_parses_and_pretty_prints_literal() {
        let expr: EmlCoordExpr = "x + ε y".parse().unwrap();

        assert_eq!(
            expr,
            EmlCoordExpr::coord(SurfaceScalarExpr::var("x"), SurfaceScalarExpr::var("y"))
        );
        assert_eq!(expr.to_string(), "x + ε y");
    }

    #[test]
    fn coord_surface_lowers_literal_to_eml_semantics() {
        let expr: EmlCoordExpr = "x + ε y".parse().unwrap();
        let lowered = expr.to_scientific();
        let expected = ScientificExpr::sub(ScientificExpr::exp(sx()), ScientificExpr::ln(sy()));

        assert_eq!(lowered, expected);

        let value = lowered
            .eval(&HashMap::from([
                ("x".to_string(), Complex::new(0.5, 0.0)),
                ("y".to_string(), Complex::new(2.0, 0.0)),
            ]))
            .unwrap();
        assert_close(value, Complex::new(0.5_f64.exp() - 2.0_f64.ln(), 0.0), 1e-12);
    }

    #[test]
    fn coord_surface_add_sub_lowers_through_scientific_expr() {
        let expr: EmlCoordExpr = "(x + ε y) - (u - ε v)".parse().unwrap();
        let lowered = expr.to_scientific();
        let expected = ScientificExpr::sub(
            ScientificExpr::sub(ScientificExpr::exp(sx()), ScientificExpr::ln(sy())),
            ScientificExpr::sub(
                ScientificExpr::exp(ScientificExpr::var("u")),
                ScientificExpr::ln(ScientificExpr::neg(ScientificExpr::var("v"))),
            ),
        );

        assert_eq!(lowered, expected);
        assert_eq!(expr.to_string(), "x + ε y - (u - ε v)");
    }

    #[test]
    fn coord_surface_scalar_multiply_is_preserved_in_ast_and_lowering() {
        let expr: EmlCoordExpr = "a * (x + ε y)".parse().unwrap();

        assert_eq!(
            expr,
            EmlCoordExpr::scalar_mul(
                SurfaceScalarExpr::var("a"),
                EmlCoordExpr::coord(SurfaceScalarExpr::var("x"), SurfaceScalarExpr::var("y"))
            )
        );
        assert_eq!(expr.to_string(), "a * (x + ε y)");

        let lowered = expr.to_scientific();
        let expected = ScientificExpr::mul(
            ScientificExpr::var("a"),
            ScientificExpr::sub(ScientificExpr::exp(sx()), ScientificExpr::ln(sy())),
        );
        assert_eq!(lowered, expected);
    }

    #[test]
    fn coord_surface_product_has_partial_support_via_lowering() {
        let expr: EmlCoordExpr = "(x + ε y) * (u + ε v)".parse().unwrap();

        assert_eq!(expr.to_string(), "(x + ε y) * (u + ε v)");
        match expr {
            EmlCoordExpr::Mul(_, _) => {}
            other => panic!("expected general multiply node, got {other:?}"),
        }

        let lowered = expr.to_scientific();
        let expected = ScientificExpr::mul(
            ScientificExpr::sub(ScientificExpr::exp(sx()), ScientificExpr::ln(sy())),
            ScientificExpr::sub(
                ScientificExpr::exp(ScientificExpr::var("u")),
                ScientificExpr::ln(ScientificExpr::var("v")),
            ),
        );
        assert_eq!(lowered, expected);
    }

    #[test]
    fn coord_surface_to_eml_uses_existing_backend() {
        let expr: EmlCoordExpr = "sin(x) + ε sqrt(y)".parse().unwrap();
        let scientific = expr.to_scientific();

        assert_eq!(expr.to_eml(), scientific.to_eml());
    }

    #[test]
    fn ball_eval_tracks_real_uncertainty() {
        let expr = ScientificExpr::add(ScientificExpr::sin(sx()), ScientificExpr::sqrt(sy()));
        let ball =
            eval_real_ball_scientific(&expr, &[("x", 0.5, 1e-6), ("y", 9.0, 1e-6)]).unwrap();

        assert!(ball.rad > 0.0);

        for x_value in [0.5 - 1e-6, 0.5, 0.5 + 1e-6] {
            for y_value in [9.0 - 1e-6, 9.0, 9.0 + 1e-6] {
                let point = eval_real_scientific(&expr, &[("x", x_value), ("y", y_value)]).unwrap();
                assert_ball_contains(ball, point);
            }
        }
    }

    #[test]
    fn ball_eval_matches_eml_tree_enclosure() {
        let expr = EmlExpr::eml(x(), y());
        let ball = eval_real_ball(&expr, &[("x", 0.2, 1e-6), ("y", 2.0, 1e-6)]).unwrap();

        for x_value in [0.2 - 1e-6, 0.2, 0.2 + 1e-6] {
            for y_value in [2.0 - 1e-6, 2.0, 2.0 + 1e-6] {
                let point = eval_real(&expr, &[("x", x_value), ("y", y_value)]).unwrap();
                assert_ball_contains(ball, point);
            }
        }
    }

    #[test]
    fn ln_ball_rejects_branch_cut_crossing() {
        let err = eval_ball_scientific(
            &ScientificExpr::ln(sx()),
            &[("x", ComplexBall::from_real(-1.0, 0.25))],
        )
        .unwrap_err();

        assert_eq!(err, BallEvalError::CrossesBranchCut("ln"));
    }

    #[test]
    fn reciprocal_ball_rejects_zero() {
        let err = eval_ball_scientific(
            &ScientificExpr::div(ScientificExpr::one(), sx()),
            &[("x", ComplexBall::from_real(0.0, 0.1))],
        )
        .unwrap_err();

        assert_eq!(err, BallEvalError::ContainsSingularity("reciprocal"));
    }
}
