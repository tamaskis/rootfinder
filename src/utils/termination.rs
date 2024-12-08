use crate::utils::convergence_data::ConvergenceData;
use crate::utils::enums::TerminationReason;
use crate::utils::solver_settings::SolverSettings;

/// Determine if the absolute step tolerance termination criterion is satisfied.
///
/// # Arguments
///
/// * `x_curr` - Current iterate.
/// * `x_next` - Next iterate.
/// * `solver_settings` - Solver settings.
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// `true` if the absolute step tolerance is satisfied, `false` otherwise.
///
/// # Note
///
/// If `solver_settings` has not defined `xatol`, then this function will just return `false`.
///
/// # Definition
///
/// The absolute step tolerance termination criterion is satisfied if
///
/// ```ignore
/// (x_next - x_curr).abs() <= xatol
/// ```
pub(crate) fn is_xatol_satisfied(
    x_curr: f64,
    x_next: f64,
    solver_settings: &SolverSettings,
    convergence_data: Option<&mut ConvergenceData>,
) -> bool {
    let satisfied = solver_settings.xatol.is_some()
        && (x_next - x_curr).abs() <= solver_settings.xatol.unwrap();
    if satisfied {
        if let Some(convergence_data) = convergence_data {
            convergence_data.termination_reason = TerminationReason::AbsoluteStepToleranceSatisfied;
        }
    }
    satisfied
}

/// Determine if the value tolerance termination criterion is satisfied.
///
/// # Arguments
///
/// * `v` - Function value.
/// * `solver_settings` - Solver settings.
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// `true` if the value tolerance is satisfied, `false` otherwise.
///
/// # Note
///
/// If `solver_settings` has not defined `vtol`, then this function will just return `false`.
///
/// # Definition
///
/// The value tolerance termination criterion is satisfied if
///
/// ```ignore
/// v.abs() <= vtol
/// ```
pub(crate) fn is_vtol_satisfied(
    v: f64,
    solver_settings: &SolverSettings,
    convergence_data: Option<&mut ConvergenceData>,
) -> bool {
    let satisfied = solver_settings.vtol.is_some() && v.abs() <= solver_settings.vtol.unwrap();
    if satisfied {
        if let Some(convergence_data) = convergence_data {
            convergence_data.termination_reason = TerminationReason::ValueToleranceSatisfied;
        }
    }
    satisfied
}

/// Determine if the absolute bracket tolerance termination criterion is satisfied.
///
/// # Arguments
///
/// * `a` - Lower bound of interval.
/// * `b` - Upper bound of interval.
/// * `solver_settings` - Solver settings.
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// `true` if the absolute bracket tolerance is satisfied, `false` otherwise.
///
/// # Note
///
/// If `solver_settings` has not defined `batol`, then this function will just return `false`.
///
/// # Note
///
/// While this crate does provide the [`crate::Interval`] struct, this function uses the lower (`a`)
/// and upper (`b`) bounds directly. This is because in the back-end of solvers, these parameters
/// are typically tracked independently as `f64`'s instead of in an [`crate::Interval`] struct.
///
/// # Definition
///
/// The absolute bracket tolerance termination criterion is satisfied if
///
/// ```ignore
/// (b - a).abs() <= batol
/// ```
#[allow(dead_code)]
pub(crate) fn is_batol_satisfied(
    a: f64,
    b: f64,
    solver_settings: &SolverSettings,
    convergence_data: Option<&mut ConvergenceData>,
) -> bool {
    let satisfied =
        solver_settings.batol.is_some() && (b - a).abs() <= solver_settings.batol.unwrap();
    if satisfied {
        if let Some(convergence_data) = convergence_data {
            convergence_data.termination_reason =
                TerminationReason::AbsoluteBracketToleranceSatisfied;
        }
    }
    satisfied
}

/// Determine if the maximum number of function evaluations termination criterion is satisfied.
///
/// # Arguments
///
/// * `n_feval` - Number of function evaluations.
/// * `solver_settings` - Solver settings.
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// `true` if the maximum number of function evaluations termination criterion is satisfied, `false`
/// otherwise.
///
/// # Note
///
/// If `solver_settings` has not defined `max_feval`, then this function will just return `false`.
///
/// # Definition
///
/// The maximum number of function evaluations termination criterion is satisfied if
///
/// ```ignore
/// n_feval >= max_feval
/// ```
pub(crate) fn is_max_feval_satisfied(
    n_feval: u32,
    solver_settings: &SolverSettings,
    convergence_data: Option<&mut ConvergenceData>,
) -> bool {
    let satisfied =
        solver_settings.max_feval.is_some() && n_feval >= solver_settings.max_feval.unwrap();
    if satisfied {
        if let Some(convergence_data) = convergence_data {
            convergence_data.termination_reason = TerminationReason::MaxFunctionEvaluationsReached;
        }
    }
    satisfied
}

#[cfg(test)]
mod is_xatol_satisfied_tests {
    use super::*;

    /// Helper function for testing [`is_xatol_satisfied`].
    ///
    /// # Arguments
    ///
    /// * `x_curr` - Current iterate.
    /// * `x_next` - Next iterate.
    /// * `xatol` - Absolute step tolerance.
    /// * `expect_satisfied` - `true` if we expect [`is_xatol_satisfied`] to return `true`,
    ///                        `false` if we expect it to return `false`.
    fn is_xatol_satisfied_test_helper(x_curr: f64, x_next: f64, xatol: f64, expect_satisied: bool) {
        // Set solver settings.
        let solver_settings = SolverSettings {
            xatol: Some(xatol),
            ..Default::default()
        };

        // Test without passing convergence data.
        let result = is_xatol_satisfied(x_curr, x_next, &solver_settings, None);
        if expect_satisied {
            assert!(result);
        } else {
            assert!(!result);
        }

        // Test with passing convergence data.
        let mut convergence_data = ConvergenceData::default();
        let result = is_xatol_satisfied(
            x_curr,
            x_next,
            &solver_settings,
            Some(&mut convergence_data),
        );
        if expect_satisied {
            assert!(result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::AbsoluteStepToleranceSatisfied
            ));
        } else {
            assert!(!result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::NotYetTerminated
            ));
        }
    }

    #[test]
    fn test_is_xatol_satisfied_default() {
        assert!(!is_xatol_satisfied(
            0.0,
            0.0,
            &SolverSettings::default(),
            None
        ));
    }

    #[test]
    fn test_is_xatol_satisfied_basic_positive_true() {
        is_xatol_satisfied_test_helper(0.5, 1.0, 0.6, true);
    }

    #[test]
    fn test_is_xatol_satisfied_basic_positive_false() {
        is_xatol_satisfied_test_helper(0.5, 1.0, 0.4, false);
    }

    #[test]
    fn test_is_xatol_satisfied_basic_negative_true() {
        is_xatol_satisfied_test_helper(1.0, 0.5, 0.6, true);
    }

    #[test]
    fn test_is_xatol_satisfied_basic_negative_false() {
        is_xatol_satisfied_test_helper(1.0, 0.5, 0.4, false);
    }

    #[test]
    fn test_is_xatol_satisfied_on_limit_positive() {
        is_xatol_satisfied_test_helper(0.5, 1.0, 0.5, true);
    }

    #[test]
    fn test_is_xatol_satisfied_on_limit_negative() {
        is_xatol_satisfied_test_helper(1.0, 0.5, 0.5, true);
    }
}

#[cfg(test)]
mod is_vtol_satisfied_tests {
    use super::*;

    /// Helper function for testing [`is_vtol_satisfied`].
    ///
    /// # Arguments
    ///
    /// * `v` - Value.
    /// * `vtol` - Value tolerance.
    /// * `expect_satisfied` - `true` if we expect [`is_vtol_satisfied`] to return `true`, `false`
    ///                        if we expect it to return `false`.
    fn is_vtol_satisfied_test_helper(v: f64, vtol: f64, expect_satisied: bool) {
        // Set solver settings.
        let solver_settings = SolverSettings {
            vtol: Some(vtol),
            ..Default::default()
        };

        // Test without passing convergence data.
        let result = is_vtol_satisfied(v, &solver_settings, None);
        if expect_satisied {
            assert!(result);
        } else {
            assert!(!result);
        }

        // Test with passing convergence data.
        let mut convergence_data = ConvergenceData::default();
        let result = is_vtol_satisfied(v, &solver_settings, Some(&mut convergence_data));
        if expect_satisied {
            assert!(result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::ValueToleranceSatisfied
            ));
        } else {
            assert!(!result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::NotYetTerminated
            ));
        }
    }

    #[test]
    fn test_is_vtol_satisfied_default() {
        assert!(!is_vtol_satisfied(0.0, &SolverSettings::default(), None));
    }

    #[test]
    fn test_is_vtol_satisfied_basic_positive_true() {
        is_vtol_satisfied_test_helper(0.5, 1.0, true);
    }

    #[test]
    fn test_is_vtol_satisfied_basic_positive_false() {
        is_vtol_satisfied_test_helper(2.0, 1.0, false);
    }

    #[test]
    fn test_is_vtol_satisfied_basic_negative_true() {
        is_vtol_satisfied_test_helper(-0.5, 1.0, true);
    }

    #[test]
    fn test_is_vtol_satisfied_basic_negative_false() {
        is_vtol_satisfied_test_helper(-2.0, 1.0, false);
    }

    #[test]
    fn test_is_vtol_satisfied_match_vtol_positive() {
        is_vtol_satisfied_test_helper(1.0, 1.0, true);
    }

    #[test]
    fn test_is_vtol_satisfied_match_vtol_negative() {
        is_vtol_satisfied_test_helper(-1.0, 1.0, true);
    }
}

#[cfg(test)]
mod is_batol_satisfied_tests {
    use super::*;

    /// Helper function for testing [`is_batol_satisfied`].
    ///
    /// # Arguments
    ///
    /// * `a` - Lower bound of interval.
    /// * `b` - Upper bound of interval.
    /// * `batol` - Absolute bracket tolerance tolerance.
    /// * `expect_satisfied` - `true` if we expect [`is_batol_satisfied`] to return `true`, `false`
    ///                        if we expect it to return `false`.
    fn is_batol_satisfied_test_helper(a: f64, b: f64, batol: f64, expect_satisied: bool) {
        // Set solver settings.
        let solver_settings = SolverSettings {
            batol: Some(batol),
            ..Default::default()
        };

        // Test without passing convergence data.
        let result = is_batol_satisfied(a, b, &solver_settings, None);
        if expect_satisied {
            assert!(result);
        } else {
            assert!(!result);
        }

        // Test with passing convergence data.
        let mut convergence_data = ConvergenceData::default();
        let result = is_batol_satisfied(a, b, &solver_settings, Some(&mut convergence_data));
        if expect_satisied {
            assert!(result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::AbsoluteBracketToleranceSatisfied
            ));
        } else {
            assert!(!result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::NotYetTerminated
            ));
        }
    }

    #[test]
    fn test_is_batol_satisfied_default() {
        assert!(!is_batol_satisfied(
            0.0,
            0.0,
            &SolverSettings::default(),
            None
        ));
    }

    #[test]
    fn test_is_batol_satisfied_basic_true() {
        is_batol_satisfied_test_helper(0.0, 1.0, 2.0, true);
    }

    #[test]
    fn test_is_batol_satisfied_basic_false() {
        is_batol_satisfied_test_helper(0.0, 1.0, 0.5, false);
    }

    #[test]
    fn test_is_batol_satisfied_match_bracket_width() {
        is_batol_satisfied_test_helper(0.0, 1.0, 1.0, true);
    }
}

#[cfg(test)]
mod is_max_feval_satisfied_tests {
    use super::*;

    /// Helper function for testing [`is_max_feval_satisfied`].
    ///
    /// # Arguments
    ///
    /// * `n_feval` - Number of function evaluations.
    /// * `max_feval` - Maximum number of function evaluations.
    /// * `expect_satisfied` - `true` if we expect [`is_max_feval_satisfied`] to return `true`,
    ///                        `false` if we expect it to return `false`.
    fn is_max_feval_satisfied_test_helper(n_feval: u32, max_feval: u32, expect_satisied: bool) {
        // Set solver settings.
        let solver_settings = SolverSettings {
            max_feval: Some(max_feval),
            ..Default::default()
        };

        // Test without passing convergence data.
        let result = is_max_feval_satisfied(n_feval, &solver_settings, None);
        if expect_satisied {
            assert!(result);
        } else {
            assert!(!result);
        }

        // Test with passing convergence data.
        let mut convergence_data = ConvergenceData::default();
        let result = is_max_feval_satisfied(n_feval, &solver_settings, Some(&mut convergence_data));
        if expect_satisied {
            assert!(result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::MaxFunctionEvaluationsReached
            ));
        } else {
            assert!(!result);
            assert!(matches!(
                convergence_data.termination_reason,
                TerminationReason::NotYetTerminated
            ));
        }
    }

    #[test]
    fn test_is_max_feval_satisfied_default() {
        assert!(!is_max_feval_satisfied(5, &SolverSettings::default(), None));
    }

    #[test]
    fn test_is_max_feval_satisfied_basic_true() {
        is_max_feval_satisfied_test_helper(4, 3, true);
    }

    #[test]
    fn test_is_max_feval_satisfied_basic_false() {
        is_max_feval_satisfied_test_helper(2, 3, false);
    }

    #[test]
    fn test_is_max_feval_satisfied_at_limit_true() {
        is_max_feval_satisfied_test_helper(3, 3, true);
    }
}
