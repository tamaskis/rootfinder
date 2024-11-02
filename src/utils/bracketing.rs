use crate::utils::convergence_data::ConvergenceData;
use crate::utils::enums::SolverError;
use crate::utils::perturb::perturb_real;
use std::fmt;

/// Interval.
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Interval {
    /// Lower bound of interval, $a\in\mathbb{R}$.
    pub a: f64,

    /// Upper bound of interval, $b\in\mathbb{R}$.
    pub b: f64,
}

impl Interval {
    /// Constructor.
    ///
    /// # Arguments
    ///
    /// `a` - Lower bound of interval, $a\in\mathbb{R}$.
    /// `b` - Upper bound of interval, $b\in\mathbb{R}$.
    ///
    /// # Note
    ///
    /// This constructor ensures that $a\leq b$.
    pub fn new(a: f64, b: f64) -> Self {
        if b > a {
            Self { a, b }
        } else {
            Self { a: b, b: a }
        }
    }

    /// Constructor from a single point.
    ///
    /// # Arguments
    ///
    /// `x` - Point, $x\in\mathbb{R}$.
    ///
    /// # Note
    ///
    /// This constructor returns the interval
    ///
    /// $[a,b]=[x,x+100\varepsilon(1+\|x\|)]$
    pub fn from_point(x: f64) -> Self {
        Self {
            a: x,
            b: perturb_real(x),
        }
    }
}

// Implementation of the std::fmt::Display trait for the Interval struct.
impl fmt::Display for Interval {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Interval(a: {}, b: {})", self.a, self.b)
    }
}

/// Find a bracketing interval.
///
/// A bracketing interval is one in which a sign change of the function occurs.
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function ($f:\mathbb{R}\to\mathbb{R}$).
/// * `ab` - Initial interval.
/// * `max_iter` - Maximum number of bracket-widening iterations. Defaults to 200.
///
/// # Returns
///
/// A result where:
///
/// * `Ok` - Bracketing interval containing a sign change in $f(x)$.
/// * `Err` - A solver error that was encountered.
pub fn bracket_sign_change(
    f: fn(f64) -> f64,
    ab: Interval,
    max_iter: Option<u32>,
    convergence_data: Option<&mut ConvergenceData>,
) -> Result<Interval, SolverError> {
    // Set maximum number of bracket-widening iterations.
    let max_iter = max_iter.unwrap_or(200);

    // If the initial interval bracket's a sign change, return it.
    if f(ab.a) * f(ab.b) < 0.0 {
        if let Some(convergence_data) = convergence_data {
            convergence_data.n_bracket_iter = 0;
            convergence_data.n_feval += 2;
        }
        return Ok(ab);
    }

    // Extract lower and upper bounds of interval and make them mutable.
    let mut a = ab.a;
    let mut b = ab.b;

    // Center of the initial interval.
    let c = (a + b) / 2.0;

    // Half-width of the initial interval.
    let mut wh = (b - a) / 2.0;

    // Variable to keep track of if sign change is found.
    let mut sign_change = false;

    // Keep expanding the interval to the left and the right until it brackets a sign change.
    let mut n_iter = 0;
    for _ in 0..max_iter {
        n_iter += 1;
        wh *= 2.0;
        a = c - wh;
        b = c + wh;
        if f(a) * f(b) < 0.0 {
            sign_change = true;
            break;
        }
    }

    // Store number of bracket-widening iterations and the number of function evaluations.
    if let Some(convergence_data) = convergence_data {
        convergence_data.n_bracket_iter = n_iter;
        convergence_data.n_feval += 2 * n_iter + 2;
    }

    // Return an error if no interval with a sign change was found.
    match sign_change {
        true => Ok(Interval::new(a, b)),
        false => Err(SolverError::BracketingIntervalNotFound),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interval_correct_order() {
        let ab = Interval::new(1.0, 2.0);
        assert_eq!(ab.a, 1.0);
        assert_eq!(ab.b, 2.0);
    }

    #[test]
    fn test_interval_incorrect_order() {
        let ab = Interval::new(2.0, 1.0);
        assert_eq!(ab.a, 1.0);
        assert_eq!(ab.b, 2.0);
    }

    #[test]
    fn test_interval_zero_width() {
        let ab = Interval::new(1.0, 1.0);
        assert_eq!(ab.a, 1.0);
        assert_eq!(ab.b, 1.0);
    }

    #[test]
    fn test_interval_print() {
        assert_eq!(
            format!("{}", Interval::new(1.0, 2.5)),
            "Interval(a: 1, b: 2.5)"
        );
    }

    #[test]
    fn test_bracket_sign_change_already_bracketed() {
        let f = |x: f64| x;
        let ab = Interval::new(-1.0, 1.0);
        let mut convergence_data = ConvergenceData::default();
        let ab_new = bracket_sign_change(f, ab, None, Some(&mut convergence_data));
        assert_eq!(ab_new.unwrap(), ab);
        assert_eq!(convergence_data.n_bracket_iter, 0);
        assert_eq!(convergence_data.n_feval, 2);
    }

    #[test]
    fn test_bracket_sign_change_from_close_initial_guess() {
        let f = |x: f64| x;
        let ab = Interval::from_point(0.0);
        let mut convergence_data = ConvergenceData::default();
        let ab_new = bracket_sign_change(f, ab, None, Some(&mut convergence_data));
        assert_eq!(
            ab_new.unwrap(),
            Interval::new(-1.1102230246251565e-14, 3.3306690738754696e-14)
        );
        assert_eq!(convergence_data.n_bracket_iter, 1);
        assert_eq!(convergence_data.n_feval, 4);
    }

    #[test]
    fn test_bracket_sign_change_from_worse_initial_guess() {
        let f = |x: f64| x;
        let ab = Interval::from_point(10.0);
        let mut convergence_data = ConvergenceData::default();
        let ab_new = bracket_sign_change(f, ab, None, Some(&mut convergence_data));
        assert_eq!(
            ab_new.unwrap(),
            Interval::new(-7.249999999999877, 27.25000000000012)
        );
        assert_eq!(convergence_data.n_bracket_iter, 47);
        assert_eq!(convergence_data.n_feval, 96);
    }

    #[test]
    fn test_bracket_sign_change_not_bracketing() {
        let f = |x: f64| x;
        let ab = Interval::new(50.0, 100.0);
        let mut convergence_data = ConvergenceData::default();
        let ab_new = bracket_sign_change(f, ab, None, Some(&mut convergence_data));
        assert_eq!(ab_new.unwrap(), Interval::new(-25.0, 175.0));
        assert_eq!(convergence_data.n_bracket_iter, 2);
        assert_eq!(convergence_data.n_feval, 6);
    }

    #[test]
    fn test_bracket_sign_change_error() {
        let f = |x: f64| x.powi(2) + 1.0;
        let ab = Interval::new(-1.0, 1.0);
        let mut convergence_data = ConvergenceData::default();
        let result = bracket_sign_change(f, ab, None, Some(&mut convergence_data));
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }
}
