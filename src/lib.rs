use std::collections::HashMap;
use std::f64::consts::{FRAC_PI_2, PI};
use std::fmt;

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
        Self::neg(Self::one())
    }

    pub fn i() -> Self {
        Self::sqrt(Self::minus_one())
    }

    pub fn pi() -> Self {
        Self::neg(Self::div(Self::ln(Self::minus_one()), Self::i()))
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

    pub fn sub(left: Self, right: Self) -> Self {
        Self::eml(Self::ln(left), Self::exp(right))
    }

    pub fn neg(expr: Self) -> Self {
        Self::sub(Self::zero(), expr)
    }

    pub fn add(left: Self, right: Self) -> Self {
        Self::sub(left, Self::neg(right))
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
        Self::exp(Self::sub(Self::ln(left), Self::ln(right)))
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
        Self::neg(Self::one())
    }

    pub fn i() -> Self {
        Self::sqrt(Self::minus_one())
    }

    pub fn pi() -> Self {
        Self::mul(Self::neg(Self::i()), Self::ln(Self::minus_one()))
    }

    pub fn sqrt(expr: Self) -> Self {
        Self::pow(expr, Self::half())
    }

    pub fn sin(expr: Self) -> Self {
        let j = Self::neg(Self::i());
        let two_i = Self::mul(Self::two(), j.clone());
        let ix = Self::mul(j, expr.clone());
        let neg_ix = Self::neg(ix.clone());
        Self::div(Self::sub(Self::exp(ix), Self::exp(neg_ix)), two_i)
    }

    pub fn cos(expr: Self) -> Self {
        let j = Self::neg(Self::i());
        let ix = Self::mul(j, expr.clone());
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
        let i = Self::i();
        let j = Self::neg(i.clone());
        let one = Self::one();
        let iz = Self::mul(j, expr.clone());
        let sqrt_term = Self::sqrt(Self::sub(one.clone(), Self::mul(expr.clone(), expr)));
        let inner = Self::add(iz, sqrt_term);
        Self::mul(i, Self::ln(inner))
    }

    pub fn acos(expr: Self) -> Self {
        let i = Self::i();
        let j = Self::neg(i.clone());
        let one = Self::one();
        let sqrt_term = Self::sqrt(Self::sub(one, Self::mul(expr.clone(), expr.clone())));
        let inner = Self::add(expr, Self::mul(j, sqrt_term));
        Self::mul(i, Self::ln(inner))
    }

    pub fn atan(expr: Self) -> Self {
        let i = Self::i();
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
