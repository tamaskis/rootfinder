use crate::utils::{
    bracketing::{bracket_sign_change, Interval},
    convergence_data::ConvergenceData,
    enums::SolverError,
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
/// The following solver settings are used by this function:
///
/// * `batol`
/// * `vtol`
/// * `max_iter`
/// * `max_feval`
/// * `rebracket`
/// * `max_bracket_iter`
/// * `store_iterates`
///
/// # Example
///
/// Finding the root of $f(x)=x^{2}-1$ in the interval $[0,\infty)$.
///
/// ```
/// use rootfinder::{root_bisection, Interval};
/// use numtest::*; // for checking result
///
/// // Define the function.
/// let f = |x: f64| x.powi(2) - 1.0;
///
/// // We want the root in the interval [0,∞). Therefore, we use an initial interval of
/// // [a,b] = [0,9999999]. Finding this root using the bisection method,
/// let result = root_bisection(f, Interval::new(0.0, 9999999.0), None, None);
/// let root = result.unwrap();
/// assert_equal_to_decimal!(root, 1.0, 16);
/// ```
pub fn root_bisection(
    f: fn(f64) -> f64,
    ab: Interval,
    solver_settings: Option<&SolverSettings>,
    mut convergence_data: Option<&mut ConvergenceData>,
) -> Result<f64, SolverError> {
    // Set solver settings.
    let solver_settings: &SolverSettings = solver_settings.unwrap_or(&DEFAULT_SOLVER_SETTINGS);

    // Default "rebracket" to "false" if not input.
    let rebracket = solver_settings.rebracket.unwrap_or(false);

    // Set the bracket tolerance.
    let batol = solver_settings.batol.unwrap_or(2.0 * f64::EPSILON);

    // New bracketing interval in case the input bracketing interval doesn't actually bracket a sign
    // change.
    let mut ab_new = Interval::new(f64::NAN, f64::NAN);

    // Rebrackets to ensure a sign change (if rebracketing is allowed).
    if rebracket {
        match bracket_sign_change(
            f,
            ab,
            solver_settings.max_bracket_iter,
            convergence_data.as_deref_mut(),
        ) {
            Ok(ab) => ab_new = ab,
            Err(e) => return Err(e),
        };
    }

    // Move the function evaluation count outside the convergence data struct (we will write it back
    // at the end).
    let mut n_feval: u32 = 0;
    if let Some(convergence_data) = convergence_data.as_deref_mut() {
        n_feval = convergence_data.n_feval;
    }

    // Get the lower and upper bounds of the initial bracketing interval.
    //  --> Note: if ab_new.a and ab_new.b aren't NaN, this is because we have re-bracketed, and we
    //            should use the "new" bracketing interval.
    let mut a;
    let mut b;
    if !ab_new.a.is_nan() && !ab_new.b.is_nan() {
        a = ab_new.a;
        b = ab_new.b;
    } else {
        a = ab.a;
        b = ab.b;
    }

    // Determine k₁⸝₂.
    let k_12 = get_k12(Interval::new(a, b), batol);

    // Set the maximum number of iterations allowed.
    let mut max_iter = solver_settings.max_iter.unwrap_or(k_12);
    max_iter = max_iter.min(k_12);

    // Initial guess.
    let mut c = (a + b) / 2.0;

    // Function evaluation at bracketing interval midpoint.
    let mut fc: f64;

    // Return the initial guess if the maximum number of function evaluations has already been
    // met/exceeded or if no iteration is required (initial interval has width within bracket
    // tolerance).
    if (solver_settings.max_feval.is_some() && n_feval >= solver_settings.max_feval.unwrap())
        || max_iter == 0
    {
        if let Some(convergence_data) = convergence_data {
            convergence_data.x_all.push(c);
            convergence_data.a_all.push(a);
            convergence_data.b_all.push(b);
            convergence_data.f_all.push(f(c));
            convergence_data.n_iter = 0;
            convergence_data.n_feval = n_feval;
        }
        return Ok(c);
    }

    // Evaluate the function at the lower bound of the bracketing interval.
    let mut fa = f(a);
    n_feval += 1;

    // Adjust the maximum number of iterations to account for the maximum number of function
    // evaluations.
    if solver_settings.max_feval.is_some() {
        max_iter = max_iter.min(solver_settings.max_feval.unwrap() - n_feval);
    }

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

    // Stores the number of function evaluations.
    if let Some(convergence_data) = convergence_data {
        convergence_data.n_feval = n_feval;
    }

    // Converged root.
    Ok(c)
}

/// Bisection method for finding the root of a univariate, scalar-valued function.
///
/// This implementation is a faster version of [`root_bisection`], with the following changes:
///
/// * It does not allow the specification of solver settings.
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
/// use rootfinder::{root_bisection_fast, Interval};
/// use numtest::*; // for checking result
///
/// // Define the function.
/// let f = |x: f64| x.powi(2) - 1.0;
///
/// // We want the root in the interval [0,∞). Therefore, we use an initial interval of
/// // [a,b] = [0,9999999]. Finding this root using the bisection method,
/// let root = root_bisection_fast(f, Interval::new(0.0, 9999999.0));
/// assert_equal_to_decimal!(root, 1.0, 16);
/// ```
pub fn root_bisection_fast(f: fn(f64) -> f64, ab: Interval) -> f64 {
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
    /// * `n_iter_exp` - Expected number of iterations.
    /// * `n_feval_exp` - Expected number of function evaluations.
    /// * `n_bracket_iter_exp` - Expected number of bracketing iterations.
    /// * `xatol` - Absolute tolerance for checking if the root (as computed by [`root_bisection`])
    ///             matches the expected root. Defaults to [`f64::EPSILON`].
    /// * `vtol` - Value tolerance for checking if the function value at the root (as computed by
    ///            [`root_bisection`]) is sufficiently close to 0. Defaults to [`f64::EPSILON`].
    #[allow(clippy::too_many_arguments)]
    fn root_bisection_test_helper(
        f: fn(f64) -> f64,
        ab: Interval,
        solver_settings: Option<&SolverSettings>,
        x_exp: f64,
        n_iter_exp: u32,
        n_feval_exp: u32,
        n_bracket_iter_exp: u32,
        xatol: Option<f64>,
        vtol: Option<f64>,
    ) {
        // Set default values.
        let xatol = xatol.unwrap_or(f64::EPSILON);
        let vtol = vtol.unwrap_or(f64::EPSILON);
        let mut convergence_data = ConvergenceData::default();

        // Solve for the root using both the "full" and "fast" versions of the bisection method.
        let x = root_bisection(f, ab, solver_settings, Some(&mut convergence_data)).unwrap();
        let x_fast = root_bisection_fast(f, ab);

        // Check results using the "full" version.
        assert_equal_to_atol!(x, x_exp, xatol);
        assert_equal_to_atol!(f(x), 0.0, vtol);

        // Check parity between full and fast implementations if not testing rebracketing (the fast
        // implementation does not rebracket).
        if solver_settings.is_some() && !solver_settings.unwrap().rebracket.unwrap_or(true) {
            assert_eq!(x, x_fast);
        }

        // Check number of iterations, function evaluations, and bracketing iterations.
        assert_eq!(convergence_data.n_iter, n_iter_exp);
        assert_eq!(convergence_data.n_feval, n_feval_exp);
        assert_eq!(convergence_data.n_bracket_iter, n_bracket_iter_exp);
    }

    #[test]
    fn test_root_bisection_root_at_midpoint() {
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 2.0),
            None,
            1.0,
            52,
            53,
            0,
            None,
            Some(2.0 * f64::EPSILON),
        );
    }

    #[test]
    fn test_root_bisection_root_at_lower_bound() {
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(1.0, 2.0),
            None,
            1.0,
            51,
            52,
            0,
            None,
            Some(2.0 * f64::EPSILON),
        );
    }

    #[test]
    fn test_root_bisection_root_at_upper_bound() {
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 1.0),
            None,
            1.0,
            51,
            52,
            0,
            None,
            Some(2.0 * f64::EPSILON),
        );
    }

    #[test]
    fn test_root_bisection_large_initial_interval() {
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 9999999.0),
            None,
            1.0,
            75,
            76,
            0,
            None,
            None,
        );
    }

    #[test]
    fn test_root_bisection_root_within_tolerance() {
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(1.0 - f64::EPSILON, 1.0 + f64::EPSILON),
            None,
            1.0,
            0,
            0,
            0,
            None,
            None,
        );
    }

    #[test]
    fn test_root_bisection_zero_bracket_width() {
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(1.0, 1.0),
            None,
            1.0,
            0,
            0,
            0,
            None,
            None,
        );
    }

    #[test]
    fn test_root_bisection_constant_function() {
        root_bisection_test_helper(
            |_x: f64| 1.0,
            Interval::new(0.0, 2.0),
            None,
            2.0,
            52,
            53,
            0,
            None,
            Some(2.0),
        );
    }

    #[test]
    fn test_root_bisection_batol() {
        let solver_settings = SolverSettings {
            batol: Some(0.1),
            ..Default::default()
        };
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            3,
            4,
            0,
            Some(solver_settings.batol.unwrap() / 2.0), // Should be within half the bracket tolerance of the true root.
            Some(0.1),
        );
    }

    #[test]
    fn test_root_bisection_vtol() {
        let solver_settings = SolverSettings {
            vtol: Some(0.01),
            ..Default::default()
        };
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            6,
            7,
            0,
            Some(solver_settings.vtol.unwrap() / 2.0),
            Some(0.01),
        );
    }

    #[test]
    fn test_root_bisection_max_iter() {
        let solver_settings = SolverSettings {
            max_iter: Some(10),
            ..Default::default()
        };
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            10,
            11,
            0,
            Some(0.0002),
            Some(0.0004),
        );
    }

    #[test]
    fn test_root_bisection_max_feval() {
        let solver_settings = SolverSettings {
            max_feval: Some(10),
            ..Default::default()
        };
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.75, 1.5),
            Some(&solver_settings),
            1.0,
            9,
            10,
            0,
            Some(0.001),
            Some(0.002),
        );
    }

    #[test]
    fn test_root_bisection_rebracket_not_needed() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(0.0, 2.0),
            Some(&solver_settings),
            1.0,
            52,
            55,
            0,
            None,
            Some(2.0 * f64::EPSILON),
        );
    }

    #[test]
    fn test_root_bisection_rebracket_successful() {
        let solver_settings = SolverSettings {
            rebracket: Some(true),
            ..Default::default()
        };
        root_bisection_test_helper(
            |x: f64| x.powi(2) - 1.0,
            Interval::new(1.5, 2.5),
            Some(&solver_settings),
            1.0,
            53,
            60,
            2,
            None,
            Some(2.0 * f64::EPSILON),
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
            |x: f64| x.powi(2) + 1.0,
            Interval::new(-2.0, 2.0),
            Some(&solver_settings),
            Some(&mut convergence_data),
        );
        assert!(result.is_err());
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
            |x: f64| x.powi(3) - 1.0,
            Interval::new(1000.0, 1100.0),
            Some(&solver_settings),
            1.0,
            63,
            76,
            5,
            Some(9.0),
            Some(99.0),
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
            |x: f64| x.powi(3) - 1.0,
            Interval::new(1000.0, 1100.0),
            Some(&solver_settings),
            None,
        );
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SolverError::BracketingIntervalNotFound
        ));
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Bisection_method
    #[test]
    fn test_root_bisection_store_iterates() {
        let f = |x: f64| x.powi(3) - x - 2.0;
        let ab = Interval::new(1.0, 2.0);
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            max_iter: Some(15),
            store_iterates: Some(true),
            ..Default::default()
        };
        let _ = root_bisection(f, ab, Some(&solver_settings), Some(&mut convergence_data)).unwrap();
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                1.5, 1.75, 1.625, 1.5625, 1.5312500, 1.5156250, 1.5234375, 1.5195313, 1.5214844,
                1.5205078, 1.5209961, 1.5212402, 1.5213623, 1.5214233, 1.5213928
            ],
            7
        );
    }
}
