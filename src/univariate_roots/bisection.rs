use crate::utils::{
    bracketing::{initial_interval_handling, Interval, IntervalResult},
    convergence_data::ConvergenceData,
    enums::{SolverError, TerminationReason},
    solver_settings::{SolverSettings, DEFAULT_SOLVER_SETTINGS},
    termination::is_vtol_satisfied,
};

/// Determine the number of bisection iterations required to converge within a bracket tolerance.
///
/// # Arguments
///
/// * `ab` - Bracketing interval.
/// * `batol` - Bracket tolerance.
///
/// # Returns
///
/// Number of bisection iterations required to converge within the bracket tolerance, `batol`.
fn get_k12(ab: Interval, batol: f64) -> u32 {
    ((ab.b - ab.a) / batol).log2().ceil() as u32
}

/// Bisection method for finding the root of a univariate, scalar-valued function.
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function, $f(x)$ ($f:\mathbb{R}\to\mathbb{R}$).
/// * `ab` - Initial bracketing interval. This interval should bracket a sign change in $f(x)$. If
///          `solver_settings.rebracket` is `true`, then this interval will be updated in the case
///          that it does not actually bracket a sign change.
/// * `solver_settings` - Solver settings. Defaults to [`DEFAULT_SOLVER_SETTINGS`].
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// A result where:
///
/// * `Ok` - Root of $f(x)$.
/// * `Err` - A solver error that was encountered.
///
/// # Note
///
/// This function always performs one more evaluation of $f(x)$ than is strictly necessary when the
/// initial interval does bracket a sign change. However, we cannot be guaranteed that the initial
/// interval does bracket a sign change, so we explicitly check for it. Therefore, instead of just
/// evaluating $f(x)$ at the lower bound ($a$) of the initial interval, this function will also
/// evaluate $f(x)$ at the upper bound ($b$) of the initial interval. Note that in addition to just
/// checking if the initial interval brackets a sign change, it also checks if a root exists at
/// either bound of the initial interval.
///
/// [`root_bisection_fast`] does not perform a function evaluation at the upper bound of the initial
/// interval. Therefore, it implicitly assumes that the initial interval brackets a sign change in
/// $f(x)$.
///
/// # Note
///
/// The following solver settings are used by this function:
///
/// * `batol`
/// * `vtol`
/// * `max_iter`
/// * `max_feval`
/// * `rebracket`
/// * `max_bracket_iter`
///
/// # Example
///
/// Finding the root of $f(x)=x^{2}-1$ in the interval $[0,\infty)$.
///
/// ```
/// use numtest::*;
///
/// use rootfinder::{root_bisection, Interval};
///
/// // Define the function.
/// let f = |x: f64| x.powi(2) - 1.0;
///
/// // We want the root in the interval [0,∞). Therefore, we use an initial interval of
/// // [a,b] = [0,9999999]. Finding this root using the bisection method,
/// let result = root_bisection(&f, Interval::new(0.0, 9999999.0), None, None);
/// let root = result.unwrap();
/// assert_equal_to_decimal!(root, 1.0, 16);
/// ```
pub fn root_bisection(
    f: &impl Fn(f64) -> f64,
    mut ab: Interval,
    solver_settings: Option<&SolverSettings>,
    mut convergence_data: Option<&mut ConvergenceData>,
) -> Result<f64, SolverError> {
    // Since we will pre-compute the maximum number of iterations based on multiple termination
    // criteria, we need to start tracking which criterion dictates the number of iterations we
    // perform.
    let mut termination_reason = TerminationReason::AbsoluteBracketToleranceSatisfied;

    // Set solver settings.
    let solver_settings: &SolverSettings = solver_settings.unwrap_or(&DEFAULT_SOLVER_SETTINGS);

    // Set the bracket tolerance.
    let batol = solver_settings.batol.unwrap_or(2.0 * f64::EPSILON);

    // Variable to store the function evaluation at the lower bound of the bracketing interval.
    let mut fa: f64;

    // Variable to track the number of evaluations of `f` performed by this function.
    let mut n_feval: u32 = 0;

    // Handling of the initial interval.
    match initial_interval_handling(f, ab, solver_settings, convergence_data.as_deref_mut()) {
        IntervalResult::Root(root) => return Ok(root),
        IntervalResult::SolverError(e) => return Err(e),
        IntervalResult::UpdatedInterval(updated_interval) => {
            ab = updated_interval.interval;
            fa = updated_interval.fa;
            n_feval += updated_interval.n_feval;
        }
    };

    // Get the lower and upper bounds of the initial bracketing interval.
    let mut a = ab.a;
    let mut b = ab.b;

    // Determine k₁⸝₂.
    let k_12 = get_k12(Interval::new(a, b), batol);

    // Set the maximum number of iterations allowed.
    let mut max_iter = solver_settings.max_iter.unwrap_or(k_12);

    // If the maximum number of iterations is greater than the number of iterations needed for
    // convergence to the absolute bracket tolerance (k_12), then update it to k_12. Otherwise, we
    // update the termination reason to reflect that the specified maximum number of iterations
    // currently dictates the solver termination.
    if max_iter >= k_12 {
        max_iter = k_12;
    } else {
        termination_reason = TerminationReason::MaxIterationsReached;
    }

    // Initial guess (i.e. initial bracketing interval midpoint).
    let mut c = (a + b) / 2.0;

    // Adjust the maximum number of iterations to account for the maximum number of function
    // evaluations.
    if solver_settings.max_feval.is_some() {
        // Number of function evaluations remaining.
        let n_feval_remaining = solver_settings.max_feval.unwrap() - n_feval;

        // If the number of function evaluations remaining is less than the maximum number of
        // iterations, we must updated the maximum number of iterations so we don't exceed the
        // maximum number of function evaluations. Additionally, we update the termination reason to
        // reflect that the maximum number of function evaluations currently dictates the solver
        // termination.
        if n_feval_remaining < max_iter {
            max_iter = n_feval_remaining;
            termination_reason = TerminationReason::MaxFunctionEvaluationsReached;
        }
    }

    // Variable to store the function evaluation at the current root estimate.
    let mut fc;

    // Iterative solution.
    for _ in 0..max_iter {
        // Evaluate the function at the current root estimate.
        fc = f(c);
        n_feval += 1;

        // Stores kth root estimate, bracketing interval, and function evaluation.
        if let Some(convergence_data) = convergence_data.as_deref_mut() {
            convergence_data.x_all.push(c);
            convergence_data.a_all.push(a);
            convergence_data.b_all.push(b);
            convergence_data.f_all.push(fc);
            convergence_data.n_iter += 1;
        }

        // Solver termination on convergence criteria.
        if is_vtol_satisfied(fc, solver_settings, convergence_data.as_deref_mut()) {
            break;
        }

        // Update the bracketing interval.
        if fa * fc > 0.0 {
            a = c;
            fa = fc;
        } else {
            b = c;
        }

        // Update the root estimate.
        c = (a + b) / 2.0;
    }

    // Stores the last root estimate, bracketing interval, the number of function evaluations, and
    // the termination reason.
    if let Some(convergence_data) = convergence_data {
        // Store data.
        convergence_data.x_all.push(c);
        convergence_data.a_all.push(a);
        convergence_data.b_all.push(b);
        convergence_data.f_all.push(f64::NAN);
        convergence_data.n_feval = n_feval;

        // Set the termination reason if it hasn't been set by one of the termination functions.
        if let TerminationReason::NotYetTerminated = convergence_data.termination_reason {
            convergence_data.termination_reason = termination_reason;
        }
    }

    // Converged root.
    Ok(c)
}

/// Bisection method for finding the root of a univariate, scalar-valued function (fast version).
///
/// This implementation is a faster version of [`root_bisection`], with the following changes:
///
/// * It does not allow the specification of solver settings.
/// * It does not check that the initial interval brackets a sign change in $f(x)$.
/// * It does not store convergence data.
/// * The iterative solver terminates when the bracketing interval is within two times the machine
///   epsilon ([`f64::EPSILON`]).
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function, $f(x)$ ($f:\mathbb{R}\to\mathbb{R}$).
/// * `ab` - Initial bracketing interval containing a sign change in $f(x)$.
///
/// # Returns
///
/// Root of $f(x)$.
///
/// # Note
///
/// This function uses a bracket tolerance of two times the machine epsilon ([`f64::EPSILON`]).
///
/// # Example
///
/// Finding the root of $f(x)=x^{2}-1$ in the interval $[0,\infty)$.
///
/// ```
/// use numtest::*;
///
/// use rootfinder::{root_bisection_fast, Interval};
///
/// // Define the function.
/// let f = |x: f64| x.powi(2) - 1.0;
///
/// // We want the root in the interval [0,∞). Therefore, we use an initial interval of
/// // [a,b] = [0,9999999]. Finding this root using the bisection method,
/// let root = root_bisection_fast(&f, Interval::new(0.0, 9999999.0));
/// assert_equal_to_decimal!(root, 1.0, 16);
/// ```
pub fn root_bisection_fast(f: &impl Fn(f64) -> f64, ab: Interval) -> f64 {
    // Determine the number of iterations needed for convergence.
    let n_iter = ((ab.b - ab.a) / (2.0 * f64::EPSILON)).log2().ceil() as u32;

    // Make a and b mutable.
    let mut a = ab.a;
    let mut b = ab.b;

    // Initial guess.
    let mut c = (a + b) / 2.0;

    // Function evaluation at current iterate.
    let mut fc;

    // Evaluate the function at the lower bound of the bracketing interval.
    let mut fa = f(a);

    // Iterative solution.
    for _ in 0..n_iter {
        // Evaluate the function at the current root estimate.
        fc = f(c);

        // Update the bracketing interval.
        if fa * fc > 0.0 {
            a = c;
            fa = fc;
        } else {
            b = c;
        }

        // Update the root estimate.
        c = (a + b) / 2.0;
    }

    // Converged root.
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::enums::TerminationReason;
    use numtest::*;

    #[test]
    fn test_get_k12() {
        assert_eq!(
            get_k12(
                Interval::new(-f64::EPSILON, f64::EPSILON),
                2.0 * f64::EPSILON
            ),
            0
        );
        assert_eq!(get_k12(Interval::new(-1.0, 1.0), 2.0 * f64::EPSILON), 52);
    }

    /// Helper function for testing [`root_bisection`].
    ///
    /// # Arguments
    ///
    /// * `f` - Univariate, scalar-valued function, $f(x)$ ($f:\mathbb{R}\to\mathbb{R}$).
    /// * `ab` - Initial bracketing interval.
    /// * `solver_settings` - Solver settings.
    /// * `x_exp` - Expected root.
    /// * `root_tol` - Absolute tolerance for checking if the root (as computed by
    ///                [`root_bisection`]) matches the expected root. Defaults to [`f64::EPSILON`].
    /// * `value_tol` - Absolute tolerance for checking if the function value at the root (as
    ///                 computed by [`root_bisection`]) is sufficiently close to 0. Defaults to
    ///                 [`f64::EPSILON`].
    /// * `n_iter_exp` - Expected number of iterations.
    /// * `n_feval_exp` - Expected number of function evaluations.
    /// * `n_bracket_iter_exp` - Expected number of bracketing iterations.
    /// * `reason_exp` - Expected termination reason.
    #[allow(clippy::too_many_arguments)]
    fn root_bisection_test_helper(
        f: &impl Fn(f64) -> f64,
        ab: Interval,
        solver_settings: Option<&SolverSettings>,
        x_exp: f64,
        root_tol: Option<f64>,
        value_tol: Option<f64>,
        n_iter_exp: u32,
        n_feval_exp: u32,
        n_bracket_iter_exp: u32,
        reason_exp: TerminationReason,
    ) {
        // Set default values.
        let root_tol = root_tol.unwrap_or(f64::EPSILON);
        let value_tol = value_tol.unwrap_or(f64::EPSILON);
        let mut convergence_data = ConvergenceData::default();

        // Solve for the root using both the "full" and "fast" versions of the bisection method.
        let x = root_bisection(f, ab, solver_settings, Some(&mut convergence_data)).unwrap();
        let x_fast = root_bisection_fast(f, ab);

        // Check results using the "full" version.
        assert_equal_to_atol!(x, x_exp, root_tol);
        assert_equal_to_atol!(f(x), 0.0, value_tol);

        // Check parity between full and fast implementations if not testing rebracketing (the fast
        // implementation does not rebracket).
        if solver_settings.is_some() && !solver_settings.unwrap().rebracket.unwrap_or(true) {
            assert_eq!(x, x_fast);
        }

        // Check number of iterations, function evaluations, and bracketing iterations.
        assert_eq!(convergence_data.n_iter, n_iter_exp);
        assert_eq!(convergence_data.n_feval, n_feval_exp);
        assert_eq!(convergence_data.n_bracket_iter, n_bracket_iter_exp);

        // Check that the termination reason matches the expected termination reason.
        assert_eq!(convergence_data.termination_reason, reason_exp);
    }

    #[test]
    fn test_root_bisection_root_at_midpoint() {
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 2.0),
            None,
            1.0,
            None,
            Some(2.0 * f64::EPSILON),
            52,
            54,
            0,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_root_at_lower_bound() {
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(1.0, 2.0),
            None,
            1.0,
            None,
            Some(2.0 * f64::EPSILON),
            0,
            2,
            0,
            TerminationReason::RootAtLowerBound,
        );
    }

    #[test]
    fn test_root_bisection_root_at_upper_bound() {
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 1.0),
            None,
            1.0,
            None,
            Some(2.0 * f64::EPSILON),
            0,
            2,
            0,
            TerminationReason::RootAtUpperBound,
        );
    }

    #[test]
    fn test_root_bisection_large_initial_interval() {
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 9999999.0),
            None,
            1.0,
            None,
            None,
            75,
            77,
            0,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_root_within_tolerance() {
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(1.0 - f64::EPSILON, 1.0 + f64::EPSILON),
            None,
            1.0,
            None,
            None,
            0,
            2,
            0,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_zero_bracket_width() {
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(1.0, 1.0),
            None,
            1.0,
            None,
            None,
            0,
            2,
            0,
            TerminationReason::RootAtLowerBound,
        );
    }

    #[test]
    fn test_root_bisection_constant_function() {
        let solver_settings = SolverSettings::default();
        let mut convergence_data = ConvergenceData::default();
        let result = root_bisection(
            &|_x: f64| 1.0,
            Interval::new(0.0, 2.0),
            Some(&solver_settings),
            Some(&mut convergence_data),
        );
        assert!(matches!(
            result.unwrap_err(),
            SolverError::IntervalDoesNotBracketSignChange
        ));
    }

    #[test]
    fn test_root_bisection_batol() {
        let solver_settings = SolverSettings {
            batol: Some(0.1),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            Some(solver_settings.batol.unwrap() / 2.0), // Should be within half the bracket tolerance of the true root.
            Some(0.1),
            3,
            5,
            0,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_vtol() {
        let solver_settings = SolverSettings {
            vtol: Some(0.01),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            Some(solver_settings.vtol.unwrap() / 2.0),
            Some(0.01),
            6,
            8,
            0,
            TerminationReason::ValueToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_max_iter() {
        let solver_settings = SolverSettings {
            max_iter: Some(10),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            Some(0.0002),
            Some(0.0004),
            10,
            12,
            0,
            TerminationReason::MaxIterationsReached,
        );
    }

    #[test]
    fn test_root_bisection_max_feval() {
        let solver_settings = SolverSettings {
            max_feval: Some(10),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            Some(0.001),
            Some(0.002),
            8,
            10,
            0,
            TerminationReason::MaxFunctionEvaluationsReached,
        );
    }

    #[test]
    fn test_root_bisection_rebracket_not_needed() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 2.0),
            Some(&solver_settings),
            1.0,
            None,
            Some(2.0 * f64::EPSILON),
            52,
            54,
            0,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_rebracket_successful() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(1.5, 2.5),
            Some(&solver_settings),
            1.0,
            None,
            Some(2.0 * f64::EPSILON),
            53,
            59,
            2,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_rebracket_failed() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        let mut convergence_data = ConvergenceData::default();
        let result = root_bisection(
            &|x: f64| x.powi(2) + 1.0,
            Interval::new(-2.0, 2.0),
            Some(&solver_settings),
            Some(&mut convergence_data),
        );
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }

    #[test]
    fn test_root_bisection_rebracket_exceeds_max_feval() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            max_feval: Some(5),
            ..Default::default()
        };
        let result = root_bisection(
            &|x: f64| x.powi(2) - 1.0,
            Interval::new(1.5, 2.5),
            Some(&solver_settings),
            None,
        );
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }

    #[test]
    fn test_root_bisection_max_bracket_iter_successful() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            max_bracket_iter: Some(5),
            ..Default::default()
        };
        root_bisection_test_helper(
            &|x: f64| x.powi(3) - 1.0,
            Interval::new(1000.0, 1100.0),
            Some(&solver_settings),
            1.0,
            Some(9.0),
            Some(99.0),
            63,
            75,
            5,
            TerminationReason::AbsoluteBracketToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_bisection_max_bracket_iter_failed() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            max_bracket_iter: Some(2),
            ..Default::default()
        };
        let result = root_bisection(
            &|x: f64| x.powi(3) - 1.0,
            Interval::new(1000.0, 1100.0),
            Some(&solver_settings),
            None,
        );
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Bisection_method
    #[test]
    fn test_root_bisection_iterates() {
        let f = |x: f64| x.powi(3) - x - 2.0;
        let ab = Interval::new(1.0, 2.0);
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            max_iter: Some(14),
            ..Default::default()
        };
        let root =
            root_bisection(&f, ab, Some(&solver_settings), Some(&mut convergence_data)).unwrap();
        assert_eq!(root, *convergence_data.x_all.last().unwrap());
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                1.5,
                1.75,
                1.625,
                1.5625,
                1.53125,
                1.515625,
                1.5234375,
                1.51953125,
                1.521484375,
                1.5205078125,
                1.52099609375,
                1.521240234375,
                1.5213623046875,
                1.52142333984375,
                1.521392822265625
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.a_all,
            [
                1.0,
                1.5,
                1.5,
                1.5,
                1.5,
                1.5,
                1.515625,
                1.515625,
                1.51953125,
                1.51953125,
                1.5205078125,
                1.52099609375,
                1.521240234375,
                1.5213623046875,
                1.5213623046875
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.b_all,
            [
                2.0,
                2.0,
                1.75,
                1.625,
                1.5625,
                1.53125,
                1.53125,
                1.5234375,
                1.5234375,
                1.521484375,
                1.521484375,
                1.521484375,
                1.521484375,
                1.521484375,
                1.52142333984375
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.f_all,
            [
                -0.125,
                1.609375,
                0.666015625,
                0.252197265625,
                0.059112548828125,
                -0.034053802490234375,
                0.012250423431396484,
                -0.010971248149871826,
                6.221756339073181e-4,
                -5.178886465728283e-3,
                -2.279443317092955e-3,
                -8.289058605441824e-4,
                -1.034331235132413e-4,
                2.593542519662151e-4,
                f64::NAN
            ],
            16
        );
    }
}
