use crate::utils::{
    convergence_data::ConvergenceData,
    enums::SolverError,
    perturb::perturb_real,
    solver_settings::SolverSettings,
    termination::{is_vtol_satisfied, is_xatol_satisfied},
};
use core::f64;
use once_cell::sync::Lazy;

/// Default Newton's method solver settings.
///
/// TODO: document if re-exporting
/// TODO: should we do this for root_bisection
pub static DEFAULT_NEWTON_SOLVER_SETTINGS: Lazy<SolverSettings> = Lazy::new(|| SolverSettings {
    max_iter: Some(200),
    xatol: Some(1e-10),
    ..Default::default()
});

/// Newton's method for finding the root of a differentiable, univariate, scalar-valued function.
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function, $f(x)$ ($f:\mathbb{R}\to\mathbb{R}$).
/// * `df` - Derivative of $f(x)$ ($f':\mathbb{R}\to\mathbb{R}$).
/// * `x0` - Initial guess for root.
/// * `solver_settings` - Solver settings. Defaults to [`DEFAULT_NEWTON_SOLVER_SETTINGS`].
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
/// The maximum number of solver iterations defaults to `200` unless otherwise provided through
/// `solver_settings`.
///
/// # Note
///
/// The following solver settings are used by this function:
///
/// * `vtol`
/// * `max_iter`
///
/// This function does not use `max_feval` or `max_deval`. This is because the number of
/// function/derivative evaluations is exactly one more than the number of iterations, so one should
/// specify the maximum number of iterations (`max_iter`) instead.
///
/// # Example
///
/// TODO
pub fn root_newton(
    f: fn(f64) -> f64,
    df: fn(f64) -> f64,
    x0: f64,
    solver_settings: Option<&SolverSettings>,
    mut convergence_data: Option<&mut ConvergenceData>,
) -> Result<f64, SolverError> {
    // Set solver settings.
    let solver_settings: &SolverSettings =
        solver_settings.unwrap_or(&DEFAULT_NEWTON_SOLVER_SETTINGS);

    // Initialize current root estimate.
    let mut x_curr = x0;

    // Declare next root estimate and function evaluation.
    let mut x_next: f64;
    let mut f_next: f64;
    let mut df_next: f64;

    // Evaluate the function and its derivative at the initial guess.
    let mut f_curr = f(x_curr);
    let mut df_curr = df(x_curr);
    if let Some(convergence_data) = convergence_data.as_deref_mut() {
        convergence_data.x_all.push(x_curr);
        convergence_data.f_all.push(f_curr);
        convergence_data.df_all.push(df_curr);
        convergence_data.n_feval += 1;
        convergence_data.n_deval += 1;
    }

    // Return the initial guess if it is a root of f(x).
    if is_vtol_satisfied(f_curr, solver_settings, convergence_data.as_deref_mut()) {
        return Ok(x_curr);
    }

    // Iterative solution.
    for _ in 0..solver_settings.max_iter.unwrap() {
        // Perturb the current root estimate if the derivative is 0, and re-evaluate the function
        // and its derivative at the perturbed root estimate.
        if df_curr == 0.0 {
            x_curr = perturb_real(x_curr);
            f_curr = f(x_curr);
            df_curr = df(x_curr);
            if let Some(convergence_data) = convergence_data.as_deref_mut() {
                convergence_data.n_feval += 1;
                convergence_data.n_deval += 1;
            }
        }

        // Update the root estimate.
        x_next = x_curr - (f_curr / df_curr);

        // Evaluate the function and its derivative at the updated root estimate.
        f_next = f(x_next);
        df_next = df(x_next);
        if let Some(convergence_data) = convergence_data.as_deref_mut() {
            convergence_data.x_all.push(x_next);
            convergence_data.f_all.push(f_next);
            convergence_data.df_all.push(df_next);
            convergence_data.n_iter += 1;
            convergence_data.n_feval += 1;
            convergence_data.n_deval += 1;
        }

        // Solver termination on convergence criteria.
        if is_xatol_satisfied(
            x_curr,
            x_next,
            solver_settings,
            convergence_data.as_deref_mut(),
        ) {
            break;
        }
        if is_vtol_satisfied(f_next, solver_settings, convergence_data.as_deref_mut()) {
            break;
        }

        // Store updated values for next iteration.
        x_curr = x_next;
        f_curr = f_next;
        df_curr = df_next;
    }

    Ok(x_curr)
}

/// Newton's method for finding the root of a differentiable, univariate, scalar-valued function
/// (fast version).
///
/// This implementation is a faster version of [`root_newton`], with the following changes:
///
/// * It only allows specification of the absolute step tolerance, and does not allow the
///   specification of general solver settings.
/// * It does not store convergence data.
///
/// # Arguments
///
/// * `f` - Univariate, scalar-valued function, $f(x)$ ($f:\mathbb{R}\to\mathbb{R}$).
/// * `df` - Derivative of $f(x)$ ($f':\mathbb{R}\to\mathbb{R}$).
/// * `x0` - Initial guess for root.
/// * `xatol` - Absolute step tolerance.    TODO: why we choose this and hard-code `max_iter`
///
/// # Returns
///
/// Root of $f(x)$.
///
/// # Note
///
/// This function uses a maximum number of iterations of 200.
///
/// # Example
///
/// TODO
pub fn root_newton_fast(f: fn(f64) -> f64, df: fn(f64) -> f64, x0: f64, xatol: Option<f64>) -> f64 {
    // Default the absolute step tolerance to `1e-10` unless otherwise specified.
    let xatol = xatol.unwrap_or(1e-10);

    // Initialize current root estimate.
    let mut x_curr = x0;

    // Declare next root estimate and function evaluation.
    let mut x_next: f64;
    let mut f_next: f64;
    let mut df_next: f64;

    // Evaluate the function and its derivative at the initial guess.
    let mut f_curr = f(x_curr);
    let mut df_curr = df(x_curr);

    // Iterative solution.
    for _ in 1..200 {
        // Perturb the current root estimate if the derivative is 0, and re-evaluate the function
        // and its derivative at the perturbed root estimate.
        if df_curr == 0.0 {
            x_curr = perturb_real(x_curr);
            f_curr = f(x_curr);
            df_curr = df(x_curr);
        }

        // Update the root estimate.
        x_next = x_curr - (f_curr / df_curr);

        // Evaluate the function and its derivative at the updated root estimate.
        f_next = f(x_next);
        df_next = df(x_next);

        // Solver termination on convergence criteria.
        if (x_next - x_curr).abs() <= xatol {
            break;
        }

        // Store updated values for next iteration.
        x_curr = x_next;
        f_curr = f_next;
        df_curr = df_next;
    }

    x_curr
}

#[cfg(test)]
mod tests {
    use super::*;
    use numtest::*;

    /// Helper function for testing [`root_newton`].
    ///
    /// # Arguments
    ///
    /// * `f` - Univariate, scalar-valued function, $f(x)$ ($f:\mathbb{R}\to\mathbb{R}$).
    /// * `df` - Derivative of $f(x)$ ($f':\mathbb{R}\to\mathbb{R}$).
    /// * `x0` - Initial guess for root.
    /// * `solver_settings` - Solver settings.
    /// * `x_exp` - Expected root.
    /// * `root_tol` - Absolute tolerance for checking if the root (as computed by
    ///                [`root_newton`]) matches the expected root. Defaults to `1e-10`.
    /// * `value_tol` - Absolute tolerance for checking if the function value at the root (as
    ///                 computed by [`root_newton`]) is sufficiently close to 0. Defaults to
    ///                 `1e-10`.
    /// * `n_iter_exp` - Expected number of iterations.
    /// * `n_feval_exp` - Expected number of function evaluations.
    /// * `n_deval_exp` - Expected number of derivative evaluations.
    #[allow(clippy::too_many_arguments)]
    fn root_newton_test_helper(
        f: fn(f64) -> f64,
        df: fn(f64) -> f64,
        x0: f64,
        solver_settings: Option<&SolverSettings>,
        x_exp: f64,
        root_tol: Option<f64>,
        value_tol: Option<f64>,
        n_iter_exp: u32,
        n_feval_exp: u32,
        n_deval_exp: u32,
    ) {
        // Set default values.
        let root_tol = root_tol.unwrap_or(1e-10);
        let value_tol = value_tol.unwrap_or(1e-10);
        let mut convergence_data = ConvergenceData::default();

        // Absolute step tolerance to use for `root_newton_fast`.
        let xatol: Option<f64> = if let Some(solver_settings) = solver_settings {
            solver_settings.xatol
        } else {
            None
        };

        // Solve for the root using both the "full" and "fast" versions of the bisection method.
        let x = root_newton(f, df, x0, solver_settings, Some(&mut convergence_data)).unwrap();
        let x_fast = root_newton_fast(f, df, x0, xatol);

        // Check results using the "full" version.
        assert_equal_to_atol!(x, x_exp, root_tol);
        assert_equal_to_atol!(f(x), 0.0, value_tol);

        // Check parity between full and fast implementations.
        assert_eq!(x, x_fast);

        // Check number of iterations, function evaluations, and bracketing iterations.
        assert_eq!(convergence_data.n_iter, n_iter_exp);
        assert_eq!(convergence_data.n_feval, n_feval_exp);
        assert_eq!(convergence_data.n_deval, n_deval_exp);
    }

    #[test]
    fn test_root_newton_near_root() {
        root_newton_test_helper(
            |x: f64| (x - 1.0).powi(3),
            |x: f64| 3.0 * (x - 1.0).powi(2),
            1.5,
            None,
            1.0,
            Some(1e-9),
            None,
            54,
            55,
            55,
        );
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Bisection_method
    #[test]
    fn test_root_newton_iterates() {
        let f = |x: f64| x.powi(3) - x.cos();
        let df = |x: f64| 3.0 * x.powi(2) + x.sin();
        let x0 = 0.5;
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            xatol: Some(1e-14),
            max_iter: Some(6),
            ..Default::default()
        };
        let _ = root_newton(
            f,
            df,
            x0,
            Some(&solver_settings),
            Some(&mut convergence_data),
        )
        .unwrap();
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                0.5,
                1.112141637097,
                0.909672693736,
                0.867263818209,
                0.865477135298,
                0.865474033111,
                0.865474033102
            ],
            12
        );
    }
}
