use crate::utils::{
    convergence_data::ConvergenceData,
    enums::{SolverError, TerminationReason},
    perturb::perturb_real,
    solver_settings::SolverSettings,
    termination::{is_vtol_satisfied, is_xatol_satisfied},
};
use core::f64;
use once_cell::sync::Lazy;

/// Default Newton's method solver settings.
///
/// | Setting | Default Value |
/// | ------- | ------------- |
/// | `max_iter` | `200` |
/// | `xatol` | `1e-10` |
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
/// Finding the root of $f(x)=x^{2}-1$ in the interval $[0,\infty)$.
///
/// ```
/// use numtest::*;
///
/// use rootfinder::root_newton;
///
/// // Define the function.
/// let f = |x: f64| x.powi(2) - 1.0;
///
/// // Define the derivative. The derivative of f(x) = x² - 1 is f'(x) = 2x.
/// let df = |x: f64| 2.0 * x;
///
/// // We want the root in the interval [0,∞). Therefore, we use an initial guess x₀ = 10. Finding
/// // this root using Newton's method,
/// let result = root_newton(&f, &df, 10.0, None, None);
/// let root = result.unwrap();
/// assert_equal_to_decimal!(root, 1.0, 11);
/// ```
pub fn root_newton(
    f: &impl Fn(f64) -> f64,
    df: &impl Fn(f64) -> f64,
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
        // Handle the edge case where the derivative is 0.
        if df_curr == 0.0 {
            // Perturb the iterate.
            x_curr = perturb_real(x_curr);

            // Re-evaluate the function and its derivative at the perturbed root estimate.
            f_curr = f(x_curr);
            df_curr = df(x_curr);
            if let Some(convergence_data) = convergence_data.as_deref_mut() {
                convergence_data.n_feval += 1;
                convergence_data.n_deval += 1;
            }

            // Terminate if the derivative is still 0 after perturbing the iterate.
            if df_curr == 0.0 {
                if let Some(convergence_data) = convergence_data.as_deref_mut() {
                    convergence_data.termination_reason = TerminationReason::ZeroDerivative;
                }
                break;
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
            x_curr = x_next;
            break;
        }
        if is_vtol_satisfied(f_next, solver_settings, convergence_data.as_deref_mut()) {
            x_curr = x_next;
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
/// * `xatol` - Absolute step tolerance.
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
/// Finding the root of $f(x)=x^{2}-1$ in the interval $[0,\infty)$.
///
/// ```
/// use numtest::*;
///
/// use rootfinder::root_newton_fast;
///
/// // Define the function.
/// let f = |x: f64| x.powi(2) - 1.0;
///
/// // Define the derivative. The derivative of f(x) = x² - 1 is f'(x) = 2x.
/// let df = |x: f64| 2.0 * x;
///
/// // We want the root in the interval [0,∞). Therefore, we use an initial guess x₀ = 10. Finding
/// // this root using Newton's method,
/// let root = root_newton_fast(&f, &df, 10.0, None);
/// assert_equal_to_decimal!(root, 1.0, 11);
/// ```
pub fn root_newton_fast(
    f: &impl Fn(f64) -> f64,
    df: &impl Fn(f64) -> f64,
    x0: f64,
    xatol: Option<f64>,
) -> f64 {
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
    for _ in 0..200 {
        // Handle the edge case where the derivative is 0.
        if df_curr == 0.0 {
            // Perturb the iterate.
            x_curr = perturb_real(x_curr);

            // Re-evaluate the function and its derivative at the perturbed root estimate.
            f_curr = f(x_curr);
            df_curr = df(x_curr);

            // Terminate if the derivative is still 0 after perturbing the iterate.
            if df_curr == 0.0 {
                break;
            }
        }

        // Update the root estimate.
        x_next = x_curr - (f_curr / df_curr);

        // Evaluate the function and its derivative at the updated root estimate.
        f_next = f(x_next);
        df_next = df(x_next);

        // Solver termination on convergence criteria.
        if (x_next - x_curr).abs() <= xatol {
            x_curr = x_next;
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
    /// * `root_tol` - Absolute tolerance for checking if the root (as computed by [`root_newton`])
    ///   matches the expected root. Defaults to `1e-10`.
    /// * `value_tol` - Absolute tolerance for checking if the function value at the root (as
    ///   computed by [`root_newton`]) is sufficiently close to 0. Defaults to `1e-10`.
    /// * `n_iter_exp` - Expected number of iterations.
    /// * `n_feval_exp` - Expected number of function evaluations.
    /// * `n_deval_exp` - Expected number of derivative evaluations.
    /// * `reason_exp` - Expected termination reason.
    #[allow(clippy::too_many_arguments)]
    fn root_newton_test_helper(
        f: &impl Fn(f64) -> f64,
        df: &impl Fn(f64) -> f64,
        x0: f64,
        solver_settings: Option<&SolverSettings>,
        x_exp: f64,
        root_tol: Option<f64>,
        value_tol: Option<f64>,
        n_iter_exp: u32,
        n_feval_exp: u32,
        n_deval_exp: u32,
        reason_exp: TerminationReason,
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

        // Check parity between full and fast implementations if custom solver settings aren't
        // specified.
        if solver_settings.is_none() {
            assert_eq!(x, x_fast);
        }

        // Check number of iterations, function evaluations, and bracketing iterations.
        assert_eq!(convergence_data.n_iter, n_iter_exp);
        assert_eq!(convergence_data.n_feval, n_feval_exp);
        assert_eq!(convergence_data.n_deval, n_deval_exp);

        // Check that the termination reason matches the expected termination reason.
        assert_eq!(convergence_data.termination_reason, reason_exp);
    }

    #[test]
    fn test_root_newton_at_root() {
        root_newton_test_helper(
            &|x: f64| (x - 1.0).powi(3),
            &|x: f64| 3.0 * (x - 1.0).powi(2),
            1.0,
            None,
            1.0,
            Some(1e-13),
            None,
            1,
            3,
            3,
            TerminationReason::AbsoluteStepToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_newton_near_root() {
        root_newton_test_helper(
            &|x: f64| (x - 1.0).powi(3),
            &|x: f64| 3.0 * (x - 1.0).powi(2),
            1.5,
            None,
            1.0,
            Some(1e-9),
            None,
            54,
            55,
            55,
            TerminationReason::AbsoluteStepToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_newton_starting_at_stationary_point() {
        root_newton_test_helper(
            &|x: f64| x.powi(2) - 1.0,
            &|x: f64| 2.0 * x,
            0.0,
            None,
            1.0,
            Some(1e-10),
            None,
            50,
            52,
            52,
            TerminationReason::AbsoluteStepToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_newton_constant_function() {
        root_newton_test_helper(
            &|_x: f64| 1.0,
            &|_x: f64| 0.0,
            1.5,
            None,
            1.0,
            Some(0.51),
            Some(1.0),
            0,
            2,
            2,
            TerminationReason::ZeroDerivative,
        );
    }

    #[test]
    fn test_root_newton_xatol() {
        let mut solver_settings = DEFAULT_NEWTON_SOLVER_SETTINGS.clone();
        solver_settings.xatol = Some(0.001);
        root_newton_test_helper(
            &|x: f64| (x - 1.0).powi(3),
            &|x: f64| 3.0 * (x - 1.0).powi(2),
            1.5,
            Some(&solver_settings),
            1.0,
            Some(0.003),
            Some(1e-7),
            14,
            15,
            15,
            TerminationReason::AbsoluteStepToleranceSatisfied,
        );
    }

    #[test]
    fn test_root_newton_vtol() {
        let mut solver_settings = DEFAULT_NEWTON_SOLVER_SETTINGS.clone();
        solver_settings.vtol = Some(0.001);
        root_newton_test_helper(
            &|x: f64| (x - 1.0).powi(3),
            &|x: f64| 3.0 * (x - 1.0).powi(2),
            1.5,
            Some(&solver_settings),
            1.0,
            Some(0.15),
            Some(0.001),
            4,
            5,
            5,
            TerminationReason::ValueToleranceSatisfied,
        );
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Newton%27s_method#Solution_of_cos(x)_=_x3_using_Newton's_method
    ///
    /// # Note
    ///
    /// Midway through, the iterates in the reference above diverge slightly from what we have here.
    #[test]
    fn test_root_newton_iterates_1() {
        let f = |x: f64| x.powi(2) - 612.0;
        let df = |x: f64| 2.0 * x;
        let x0 = 1.0;
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            xatol: Some(1e-14),
            max_iter: Some(10),
            ..Default::default()
        };
        let root = root_newton(
            &f,
            &df,
            x0,
            Some(&solver_settings),
            Some(&mut convergence_data),
        )
        .unwrap();
        assert_eq!(root, *convergence_data.x_all.last().unwrap());
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                1.0,
                306.5,
                154.2483686786297,
                79.10799786435472,
                43.42212868215148,
                28.758162428779126,
                25.019538536995714,
                24.74021067122501,
                24.738633803961573,
                24.738633753705965,
                24.73863375370596
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.f_all,
            [
                -6.11e2,
                93330.25,
                23180.559240018472,
                5646.07532610675,
                1273.4812592893222,
                215.03190628004336,
                13.977308604213704,
                0.07802405659595024,
                2.486510197741154e-6,
                1.1368683772161603e-13,
                -1.1368683772161603e-13
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.df_all,
            [
                2.0,
                613.0,
                308.4967373572594,
                158.21599572870943,
                86.84425736430296,
                57.51632485755825,
                50.03907707399143,
                49.48042134245002,
                49.477267607923146,
                49.47726750741193,
                49.47726750741192
            ],
            16
        );
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Newton%27s_method#Solution_of_cos(x)_=_x3_using_Newton's_method
    #[test]
    fn test_root_newton_iterates_2() {
        let f = |x: f64| x.powi(2) - 612.0;
        let df = |x: f64| 2.0 * x;
        let x0 = 10.0;
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            xatol: Some(1e-14),
            max_iter: Some(5),
            ..Default::default()
        };
        let root = root_newton(
            &f,
            &df,
            x0,
            Some(&solver_settings),
            Some(&mut convergence_data),
        )
        .unwrap();
        assert_eq!(root, *convergence_data.x_all.last().unwrap());
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                10.0,
                35.6,
                26.395505617977527,
                24.790635492455475,
                24.738688294075324,
                24.738633753766084
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.f_all,
            [
                -512.0,
                655.3600000000001,
                84.72271682868313,
                2.5756081197931735,
                2.698511419453098e-3,
                2.9746161089860834e-9
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.df_all,
            [
                20.0,
                71.2,
                52.791011235955054,
                49.58127098491095,
                49.47737658815065,
                49.47726750753217
            ],
            16
        );
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Newton%27s_method#Solution_of_cos(x)_=_x3_using_Newton's_method
    ///
    /// # Warning
    ///
    /// In the reference above, they have f(x₄) = 6.1424 ✕ 10⁻¹³. However, if you use the unrounded
    /// root estimate, you get f(x₄) = 6.8212 ✕ 10⁻¹³.
    #[test]
    fn test_root_newton_iterates_3() {
        let f = |x: f64| x.powi(2) - 612.0;
        let df = |x: f64| 2.0 * x;
        let x0 = -20.0;
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            xatol: Some(1e-14),
            max_iter: Some(4),
            ..Default::default()
        };
        let root = root_newton(
            &f,
            &df,
            x0,
            Some(&solver_settings),
            Some(&mut convergence_data),
        )
        .unwrap();
        assert_eq!(root, *convergence_data.x_all.last().unwrap());
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                -20.0,
                -25.3,
                -24.744861660079053,
                -24.738634537440753,
                -24.738633753705976
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.f_all,
            [
                -212.0,
                28.090000000000032,
                0.3081785764502456,
                3.877705648847041e-5,
                6.821210263296962e-13
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.df_all,
            [
                -40.0,
                -50.6,
                -49.489723320158106,
                -49.477269074881505,
                -49.47726750741195
            ],
            16
        );
    }

    /// # References
    ///
    /// * https://en.wikipedia.org/wiki/Newton%27s_method#Solution_of_cos(x)_=_x3_using_Newton's_method
    #[test]
    fn test_root_newton_iterates_4() {
        let f = |x: f64| x.powi(3) - x.cos();
        let df = |x: f64| 3.0 * x.powi(2) + x.sin();
        let x0 = 0.5;
        let mut convergence_data = ConvergenceData::default();
        let solver_settings = SolverSettings {
            xatol: Some(1e-14),
            max_iter: Some(6),
            ..Default::default()
        };
        let root = root_newton(
            &f,
            &df,
            x0,
            Some(&solver_settings),
            Some(&mut convergence_data),
        )
        .unwrap();
        assert_eq!(root, *convergence_data.x_all.last().unwrap());
        assert_arrays_equal_to_decimal!(
            convergence_data.x_all,
            [
                0.5,
                1.1121416370972725,
                0.9096726937368068,
                0.8672638182088165,
                0.8654771352982646,
                0.8654740331109566,
                0.8654740331016144
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.f_all,
            [
                -0.7525825618903728,
                0.9328201795040982,
                0.13875403935061037,
                0.005393998041341108,
                9.333106352094056e-6,
                2.8106295069108e-11,
                -2.220446049250313e-16
            ],
            16
        );
        assert_arrays_equal_to_decimal!(
            convergence_data.df_all,
            [
                1.229425538604203,
                4.6072259973390155,
                3.2718160437673296,
                3.019001306546996,
                3.0085566812522133e0,
                3.008538560967953,
                3.008538560913384
            ],
            16
        );
    }
}
