use crate::utils::convergence_data::ConvergenceData;
use crate::utils::enums::SolverError;
use crate::utils::enums::TerminationReason;
use crate::utils::perturb::perturb_real;
use crate::utils::solver_settings::SolverSettings;
use crate::utils::termination::{is_btol_satisfied, is_vtol_satisfied};
use std::fmt;

/// Updated interval with some associated metadata for bracketing methods.
///
/// # Note
///
/// The primary purpose of this struct is to store data produced by [`initial_interval_handling`].
#[derive(Debug, PartialEq)]
pub(crate) struct UpdatedInterval {
    /// Updated interval.
    pub(crate) interval: Interval,

    /// Function evaluation at the lower bound of the updated interval.
    pub(crate) fa: f64,

    /// Number of function evaluations performed to get the updated interval.
    pub(crate) n_feval: u32,
}

impl UpdatedInterval {
    /// Constructor.
    ///
    /// # Arguments
    ///
    /// * `interval` - Updated interval.
    /// * `fa` - Function evaluation at the lower bound of the updated interval.
    /// * `n_feval` - Number of function evaluations performed to get the updated interval.
    ///
    /// # Returns
    ///
    /// Updated interval with some associated metadata for bracketing methods.
    pub(crate) fn new(interval: Interval, fa: f64, n_feval: u32) -> UpdatedInterval {
        UpdatedInterval {
            interval,
            fa,
            n_feval,
        }
    }
}

#[derive(Debug)]
pub enum IntervalResult {
    /// The updated interval produced by [`initial_interval_handling`] (if a root wasn't found or a
    /// solver error wasn't encountered).
    UpdatedInterval(UpdatedInterval),

    /// The root of the function if one was found during the handling of the initial interval.
    Root(f64),

    /// A solver error that was encountered during the handling of the initial interval.
    ///
    /// # Note
    ///
    /// This can be either [`SolverError::BracketingIntervalNotFound`] if rebracketing was attempted
    /// but failed, or [`SolverError::IntervalDoesNotBracketSignChange`] if the initial interval did
    /// not bracket a sign change and no rebracketing was attempted.
    SolverError(SolverError),
}

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
    /// This constructor ensures that $a<b$.
    ///
    /// * If $a=b$, then $b$ is set to $b+100\varepsilon(1+\|b\|)$ (using `perturb_real`).
    /// * If $a>b$, then $a$ and $b$ are swapped.
    pub fn new(a: f64, b: f64) -> Self {
        if a == b {
            Self {
                a,
                b: perturb_real(a),
            }
        } else if b > a {
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
    /// $\[a,b\]=[x,x+100\varepsilon(1+\|x\|)]$
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
/// A bracketing interval is one in which a sign change of the function, $f(x)$, occurs.
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function ($f:\mathbb{R}\to\mathbb{R}$).
/// * `ab` - Initial interval.
/// * `fa` - Function evaluation at the lower bound of the initial interval, $f(a)$.
/// * `fb` - Function evaluation at the upper bound of the initial interval, $f(b)$.
/// * `max_bracket_iter` - Maximum number of iterations to find a bracketing interval allowed.
///   Defaults to 200.
///
/// # Returns
///
/// A result where:
///
/// * `Ok` - A tuple where:
///
///     * `0` - Bracketing interval containing a sign change in $f(x)$.
///     * `1` - Function evaluation at the lower bound of the returned bracketing interval, $f(a)$.
///     * `2` - Number of iterations to find a bracketing interval.
///     * `3` - Number of evaluations of `f` performed by this function.
///
/// * `Err` - A solver error that was encountered.
///
/// # Note
///
/// Each bracket
pub fn bracket_sign_change(
    f: &impl Fn(f64) -> f64,
    ab: Interval,
    mut fa: f64,
    mut fb: f64,
    max_bracket_iter: u32,
) -> Result<(Interval, f64, u32, u32), SolverError> {
    // If the initial interval bracket's a sign change, return it.
    if fa * fb < 0.0 {
        return Ok((ab, fa, 0, 0));
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

    // Variables to track the number of iterations to find a bracketing interval and the number of
    // evaluations of `f` performed by this function.
    let mut n_bracket_iter: u32 = 0;
    let mut n_feval: u32 = 0;

    // Keep expanding the interval to the left and the right until it brackets a sign change.
    for _ in 0..max_bracket_iter {
        // Track the number of iterations.
        n_bracket_iter += 1;

        // Increase the interval half-width.
        wh *= 2.0;
        a = c - wh;
        b = c + wh;

        // Evaluate the function at the new interval bounds.
        fa = f(a);
        fb = f(b);
        n_feval += 2;

        // Determine if the new interval brackets a sign change.
        if fa * fb < 0.0 {
            sign_change = true;
            break;
        }
    }

    // Return an error if no interval with a sign change was found.
    match sign_change {
        true => Ok((Interval::new(a, b), fa, n_bracket_iter, n_feval)),
        false => Err(SolverError::BracketingIntervalNotFound),
    }
}

/// Handling of the initial interval passed to a bracketing method.
///
/// This function performs the following checks on the initial interval in the presented order:
///
/// 1. Checks for a root at the lower bound of the initial interval, using `vtol` if specified in
///    `solver_settings`.
/// 2. Checks for a root at the upper bound of the initial interval, using `vtol` if specified in
///    `solver_settings`.
/// 3. If the initial interval brackets a sign change, checks if the absolute bracket tolerance is
///    satisfied, in which case the midpoint of the initial interval is a root.
/// 4. If the initial interval brackets a sign change, checks if the relative bracket tolerance is
///    satisfied, in which case the midpoint of the initial interval is a root.
///
/// If none of the checks above yield a root or an error, this function then tries to find a
/// bracketing interval (i.e. an interval in which a sign change of the function, $f(x)$, occurs.)
/// if `solver_settings.rebracket` is `true`.
///
/// Finally, this function returns the updated interval (with some associated metadata), or an error
/// if no rebracketing was done and the initial interval does not bracket a sign change.
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function ($f:\mathbb{R}\to\mathbb{R}$).
/// * `ab` - Initial interval.
/// * `solver_settings` - Solver settings.
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// The result of the initial interval handling. See [`IntervalResult`] for more detail.
///
/// # Warning
///
/// No evaluations of `f` should have been performed before calling this function.
///
/// # Note
///
/// If `solver_settings.max_feval` is specified, then the maximum number of iterations allowed to
/// find a bracketing interval is adjusted accordingly.
///
/// # Note
///
/// By default:
///
/// * This function does not perform any rebracketing unless otherwise specified through
///   `solver_settings`.
/// * This function allows `bracket_sign_change` to perform up to and including 200 iterations to
///   find a bracketing interval if `solver_settings.rebracket` is `true`, unless some other limit
///   is specified by `solver_settings.max_bracket_iter`.
pub fn initial_interval_handling(
    f: &impl Fn(f64) -> f64,
    ab: Interval,
    solver_settings: &SolverSettings,
    mut convergence_data: Option<&mut ConvergenceData>,
) -> IntervalResult {
    // Variable to track the number of evaluations of `f` performed by this function.
    let mut n_feval: u32 = 0;

    // Default `rebracket` to `false` unless otherwise specified by `solver_settings`.
    let rebracket = solver_settings.rebracket.unwrap_or(false);

    // Evaluate the function at the bounds of the initial interval.
    let mut fa = f(ab.a);
    let fb = f(ab.b);
    n_feval += 2;
    if let Some(convergence_data) = convergence_data.as_deref_mut() {
        convergence_data.n_feval += n_feval;
    }

    // Determine if there is a root at either bound of the interval.
    let root_at_lower_bound =
        is_vtol_satisfied(fa, solver_settings, convergence_data.as_deref_mut()) || fa == 0.0;
    let root_at_upper_bound =
        is_vtol_satisfied(fb, solver_settings, convergence_data.as_deref_mut()) || fb == 0.0;

    // Handling for the case where there is a root at one of the bounds.
    if root_at_lower_bound || root_at_upper_bound {
        // Set the root, prioritizing the root at the lower bounds.
        let root = if root_at_lower_bound { ab.a } else { ab.b };

        // Store the convergence data.
        if let Some(convergence_data) = convergence_data.as_deref_mut() {
            convergence_data.x_all.push(root);
            convergence_data.a_all.push(ab.a);
            convergence_data.b_all.push(ab.b);
            if root_at_lower_bound {
                convergence_data.f_all.push(fa);
                if convergence_data.termination_reason != TerminationReason::ValueToleranceSatisfied
                {
                    convergence_data.termination_reason = TerminationReason::RootAtLowerBound;
                }
            } else {
                convergence_data.f_all.push(fb);
                if convergence_data.termination_reason != TerminationReason::ValueToleranceSatisfied
                {
                    convergence_data.termination_reason = TerminationReason::RootAtUpperBound;
                }
            }
        }

        // Return the root.
        return IntervalResult::Root(root);
    }

    // If the interval brackets a sign change and the bracket tolerance(s) are satisfied, then the
    // midpoint of the interval is the root.
    if (fa.signum() != fb.signum())
        && is_btol_satisfied(ab.a, ab.b, solver_settings, convergence_data.as_deref_mut())
    {
        let root = (ab.a + ab.b) / 2.0;
        if let Some(convergence_data) = convergence_data.as_deref_mut() {
            convergence_data.x_all.push(root);
            convergence_data.a_all.push(ab.a);
            convergence_data.b_all.push(ab.b);
            convergence_data.f_all.push(f64::NAN);
        }
        return IntervalResult::Root(root);
    }

    // Perform rebracketing if requested.
    if rebracket {
        // Set the maximum number of iterations allowed for finding a bracketing interval.
        let mut max_bracket_iter = solver_settings.max_bracket_iter.unwrap_or(200);

        // Determine the number of function evaluations remaining.
        let n_feval_remaining = solver_settings
            .max_feval
            .map(|max_feval| max_feval - n_feval);

        // Update the number of iterations allowed to find a bracketing interval to account for the
        // number of remaining function evaluations.
        if let Some(n_feval_remaining) = n_feval_remaining {
            max_bracket_iter = max_bracket_iter.min(n_feval_remaining / 2);
        }

        // Find a bracketing interval, or return the appropriate solver error.
        match bracket_sign_change(f, ab, fa, fb, max_bracket_iter) {
            Ok((ab, fa_, n_bracket_iter, n_feval_rebracket)) => {
                fa = fa_;
                n_feval += n_feval_rebracket;
                if let Some(convergence_data) = convergence_data {
                    convergence_data.n_bracket_iter = n_bracket_iter;
                    convergence_data.n_feval += n_feval_rebracket;
                }
                return IntervalResult::UpdatedInterval(UpdatedInterval::new(ab, fa, n_feval));
            }
            Err(err) => return IntervalResult::SolverError(err),
        }
    }
    // If not rebracketing, still check whether the interval is a bracketing interval by checking
    // the signs of f(a) and f(b).
    if fa.signum() == fb.signum() {
        IntervalResult::SolverError(SolverError::IntervalDoesNotBracketSignChange)
    } else {
        IntervalResult::UpdatedInterval(UpdatedInterval::new(ab, fa, n_feval))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::solver_settings::DEFAULT_SOLVER_SETTINGS;
    use numtest::*;

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
        assert_eq!(ab.b, 1.0000000000000444);
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
        let fa = f(ab.a);
        let fb = f(ab.b);
        let (ab_new, fa, n_bracket_iter, n_feval) =
            bracket_sign_change(&f, ab, fa, fb, 200).unwrap();
        assert_eq!(ab_new, ab);
        assert_eq!(fa, -1.0);
        assert_eq!(n_bracket_iter, 0);
        assert_eq!(n_feval, 0);
    }

    #[test]
    fn test_bracket_sign_change_from_close_initial_guess() {
        let f = |x: f64| x;
        let ab = Interval::from_point(0.0);
        let fa = f(ab.a);
        let fb = f(ab.b);
        let (ab_new, fa, n_bracket_iter, n_feval) =
            bracket_sign_change(&f, ab, fa, fb, 200).unwrap();
        assert_eq!(
            ab_new,
            Interval::new(-1.1102230246251565e-14, 3.3306690738754696e-14)
        );
        assert_eq!(fa, -1.1102230246251565e-14);
        assert_eq!(n_bracket_iter, 1);
        assert_eq!(n_feval, 2);
    }

    #[test]
    fn test_bracket_sign_change_from_worse_initial_guess() {
        let f = |x: f64| x;
        let ab = Interval::from_point(10.0);
        let fa = f(ab.a);
        let fb = f(ab.b);
        let (ab_new, fa, n_bracket_iter, n_feval) =
            bracket_sign_change(&f, ab, fa, fb, 200).unwrap();
        assert_eq!(ab_new, Interval::new(-7.249999999999877, 27.25000000000012));
        assert_eq!(fa, -7.249999999999877);
        assert_eq!(n_bracket_iter, 47);
        assert_eq!(n_feval, 94);
    }

    #[test]
    fn test_bracket_sign_change_not_bracketing() {
        let f = |x: f64| x;
        let ab = Interval::new(50.0, 100.0);
        let fa = f(ab.a);
        let fb = f(ab.b);
        let (ab_new, fa, n_bracket_iter, n_feval) =
            bracket_sign_change(&f, ab, fa, fb, 200).unwrap();
        assert_eq!(ab_new, Interval::new(-25.0, 175.0));
        assert_eq!(fa, -25.0);
        assert_eq!(n_bracket_iter, 2);
        assert_eq!(n_feval, 4);
    }

    #[test]
    fn test_bracket_sign_change_no_interval_found() {
        let f = |x: f64| x.powi(2) + 1.0;
        let ab = Interval::new(-1.0, 1.0);
        let fa = f(ab.a);
        let fb = f(ab.b);
        let result = bracket_sign_change(&f, ab, fa, fb, 200);
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }

    #[test]
    fn test_bracket_sign_change_max_bracket_iter_reached() {
        let f = |x: f64| x;
        let ab = Interval::new(10.0, 10.1);
        let fa = f(ab.a);
        let fb = f(ab.b);
        let result = bracket_sign_change(&f, ab, fa, fb, 7);
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }

    #[test]
    fn test_initial_interval_handling_already_bracketed() {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(-1.0, 1.0);

        // Check once with rebracketing and once without.
        let solver_settings_1 = SolverSettings::default();
        let solver_settings_2 = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        for solver_settings in [solver_settings_1, solver_settings_2] {
            let mut convergence_data = ConvergenceData::default();
            let interval_result =
                initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
            match interval_result {
                IntervalResult::UpdatedInterval(updated_interval) => {
                    assert_eq!(
                        updated_interval,
                        UpdatedInterval::new(Interval::new(-1.0, 1.0), -1.0, 2)
                    );
                    assert_eq!(convergence_data.x_all, vec![]);
                    assert_eq!(convergence_data.a_all, vec![]);
                    assert_eq!(convergence_data.b_all, vec![]);
                    assert_eq!(convergence_data.f_all, vec![]);
                    assert_eq!(convergence_data.n_feval, 2);
                    assert_eq!(convergence_data.n_bracket_iter, 0);
                }
                _ => panic!("Test failed."),
            }
        }
    }

    #[test]
    fn test_initial_interval_handling_rebracketing() {
        let f = |x: f64| x;
        let ab = Interval::new(10.0, 10.1);
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        let mut convergence_data = ConvergenceData::default();
        let interval_result =
            initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
        match interval_result {
            IntervalResult::UpdatedInterval(updated_interval) => {
                assert_eq!(
                    updated_interval,
                    UpdatedInterval::new(
                        Interval::new(-2.749999999999954, 22.849999999999955),
                        -2.749999999999954,
                        18
                    )
                );
                assert_eq!(convergence_data.x_all, vec![]);
                assert_eq!(convergence_data.a_all, vec![]);
                assert_eq!(convergence_data.b_all, vec![]);
                assert_eq!(convergence_data.f_all, vec![]);
                assert_eq!(convergence_data.n_feval, 18);
                assert_eq!(convergence_data.n_bracket_iter, 8);
            }
            _ => panic!("Test failed."),
        }
    }

    #[test]
    fn test_initial_interval_handling_root_at_lower_bound_no_vtol() {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(0.0, 1.0);

        // Check once with rebracketing and once without.
        let solver_settings_1 = SolverSettings::default();
        let solver_settings_2 = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        for solver_settings in [solver_settings_1, solver_settings_2] {
            let mut convergence_data = ConvergenceData::default();
            let interval_result =
                initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
            match interval_result {
                IntervalResult::Root(root) => {
                    assert_eq!(root, 0.0);
                    assert_eq!(convergence_data.x_all, vec![0.0]);
                    assert_eq!(convergence_data.a_all, vec![0.0]);
                    assert_eq!(convergence_data.b_all, vec![1.0]);
                    assert_eq!(convergence_data.f_all, vec![0.0]);
                    assert_eq!(convergence_data.n_feval, 2);
                    assert_eq!(convergence_data.n_bracket_iter, 0);
                    assert_eq!(
                        convergence_data.termination_reason,
                        TerminationReason::RootAtLowerBound
                    );
                }
                _ => panic!("Test failed."),
            }
        }
    }

    #[test]
    fn test_initial_interval_handling_root_at_lower_bound_with_vtol() {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(0.1, 1.0);

        // Check once with rebracketing and once without.
        let solver_settings_1 = SolverSettings {
            vtol: Some(0.1),
            ..Default::default()
        };
        let solver_settings_2 = SolverSettings {
            vtol: Some(0.1),
            rebracket: Some(true),
            ..Default::default()
        };
        for solver_settings in [solver_settings_1, solver_settings_2] {
            let mut convergence_data = ConvergenceData::default();
            let interval_result =
                initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
            match interval_result {
                IntervalResult::Root(root) => {
                    assert_eq!(root, 0.1);
                    assert_eq!(convergence_data.x_all, vec![0.1]);
                    assert_eq!(convergence_data.a_all, vec![0.1]);
                    assert_eq!(convergence_data.b_all, vec![1.0]);
                    assert_eq!(convergence_data.f_all, vec![0.1]);
                    assert_eq!(convergence_data.n_feval, 2);
                    assert_eq!(convergence_data.n_bracket_iter, 0);
                    assert_eq!(
                        convergence_data.termination_reason,
                        TerminationReason::ValueToleranceSatisfied
                    );
                }
                _ => panic!("Test failed."),
            }
        }
    }

    #[test]
    fn test_initial_interval_handling_root_at_upper_bound_no_vtol() {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(-1.0, 0.0);

        // Check once with rebracketing and once without.
        let solver_settings_1 = SolverSettings::default();
        let solver_settings_2 = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        for solver_settings in [solver_settings_1, solver_settings_2] {
            let mut convergence_data = ConvergenceData::default();
            let interval_result =
                initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
            match interval_result {
                IntervalResult::Root(root) => {
                    assert_eq!(root, 0.0);
                    assert_eq!(convergence_data.x_all, vec![0.0]);
                    assert_eq!(convergence_data.a_all, vec![-1.0]);
                    assert_eq!(convergence_data.b_all, vec![0.0]);
                    assert_eq!(convergence_data.f_all, vec![0.0]);
                    assert_eq!(convergence_data.n_feval, 2);
                    assert_eq!(convergence_data.n_bracket_iter, 0);
                    assert_eq!(
                        convergence_data.termination_reason,
                        TerminationReason::RootAtUpperBound
                    );
                }
                _ => panic!("Test failed."),
            }
        }
    }

    #[test]
    fn test_initial_interval_handling_root_at_upper_bound_with_vtol() {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(-1.0, -0.1);

        // Check once with rebracketing and once without.
        let solver_settings_1 = SolverSettings {
            vtol: Some(0.1),
            ..Default::default()
        };
        let solver_settings_2 = SolverSettings {
            vtol: Some(0.1),
            rebracket: Some(true),
            ..Default::default()
        };
        for solver_settings in [solver_settings_1, solver_settings_2] {
            let mut convergence_data = ConvergenceData::default();
            let interval_result =
                initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
            match interval_result {
                IntervalResult::Root(root) => {
                    assert_eq!(root, -0.1);
                    assert_eq!(convergence_data.x_all, vec![-0.1]);
                    assert_eq!(convergence_data.a_all, vec![-1.0]);
                    assert_eq!(convergence_data.b_all, vec![-0.1]);
                    assert_eq!(convergence_data.f_all, vec![-0.1]);
                    assert_eq!(convergence_data.n_feval, 2);
                    assert_eq!(convergence_data.n_bracket_iter, 0);
                    assert_eq!(
                        convergence_data.termination_reason,
                        TerminationReason::ValueToleranceSatisfied
                    );
                }
                _ => panic!("Test failed."),
            }
        }
    }

    #[test]
    fn test_initial_interval_handling_root_batol_satisfied_with_sign_change() {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(-0.1, 0.1);

        // Check once with rebracketing and once without.
        let solver_settings_1 = SolverSettings {
            batol: Some(0.2),
            ..Default::default()
        };
        let solver_settings_2 = SolverSettings {
            batol: Some(0.2),
            rebracket: Some(true),
            ..Default::default()
        };
        for solver_settings in [solver_settings_1, solver_settings_2] {
            let mut convergence_data = ConvergenceData::default();
            let interval_result =
                initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
            match interval_result {
                IntervalResult::Root(root) => {
                    assert_eq!(root, 0.0);
                    assert_eq!(convergence_data.x_all, vec![0.0]);
                    assert_eq!(convergence_data.a_all, vec![-0.1]);
                    assert_eq!(convergence_data.b_all, vec![0.1]);
                    assert_arrays_equal!(convergence_data.f_all, [f64::NAN]);
                    assert_eq!(convergence_data.n_feval, 2);
                    assert_eq!(convergence_data.n_bracket_iter, 0);
                    assert_eq!(
                        convergence_data.termination_reason,
                        TerminationReason::AbsoluteBracketToleranceSatisfied
                    );
                }
                _ => panic!("Test failed."),
            }
        }
    }

    #[test]
    fn test_initial_interval_handling_root_batol_satisfied_without_sign_change_without_rebracketing()
     {
        // Common inputs for both checks in this test.
        let f = |x: f64| x;
        let ab = Interval::new(0.1, 0.2);

        // Check once with rebracketing and once without.
        let solver_settings = SolverSettings {
            batol: Some(0.2),
            ..Default::default()
        };
        let mut convergence_data = ConvergenceData::default();
        let interval_result =
            initial_interval_handling(&f, ab, &solver_settings, Some(&mut convergence_data));
        match interval_result {
            IntervalResult::SolverError(err) => {
                assert!(matches!(err, SolverError::IntervalDoesNotBracketSignChange));
            }
            _ => panic!("Test failed."),
        }
    }

    #[test]
    fn test_initial_interval_handling_not_rebracketing_no_sign_change() {
        let f = |x: f64| x;
        let ab = Interval::new(50.0, 100.0);
        let solver_settings = &DEFAULT_SOLVER_SETTINGS;
        let interval_result = initial_interval_handling(&f, ab, solver_settings, None);
        match interval_result {
            IntervalResult::SolverError(err) => {
                assert!(matches!(err, SolverError::IntervalDoesNotBracketSignChange));
            }
            _ => panic!("Test failed."),
        }
    }

    #[test]
    fn test_initial_interval_handling_rebracketing_no_interval_found() {
        let f = |x: f64| x.powi(2) + 1.0;
        let ab = Interval::new(-1.0, 1.0);
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        let interval_result = initial_interval_handling(&f, ab, &solver_settings, None);
        match interval_result {
            IntervalResult::SolverError(err) => {
                assert!(matches!(err, SolverError::BracketingIntervalNotFound));
            }
            _ => panic!("Test failed."),
        }
    }

    #[test]
    fn test_initial_interval_handling_even_max_feval_reached_during_rebracketing() {
        let f = |x: f64| x;
        let ab = Interval::new(10.0, 10.1);
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            max_feval: Some(16),
            ..Default::default()
        };
        let interval_result = initial_interval_handling(&f, ab, &solver_settings, None);
        match interval_result {
            IntervalResult::SolverError(err) => {
                assert!(matches!(err, SolverError::BracketingIntervalNotFound));
            }
            _ => panic!("Test failed."),
        }
    }

    #[test]
    fn test_initial_interval_handling_odd_max_feval_reached_during_rebracketing() {
        let f = |x: f64| x;
        let ab = Interval::new(10.0, 10.1);
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            max_feval: Some(17),
            ..Default::default()
        };
        let interval_result = initial_interval_handling(&f, ab, &solver_settings, None);
        match interval_result {
            IntervalResult::SolverError(err) => {
                assert!(matches!(err, SolverError::BracketingIntervalNotFound));
            }
            _ => panic!("Test failed."),
        }
    }
}
