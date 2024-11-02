use crate::utils::convergence_data::ConvergenceData;
use crate::utils::enums::TerminationReason;
use crate::utils::solver_settings::SolverSettings;

/// Determine if the absolute bracket tolerance is satisfied.
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
/// * `true` if the absolute bracket tolerance is satisfied, `false` otherwise.
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

/// Determine if the value tolerance is satisfied.
///
/// # Arguments
///
/// * `v` - Function value.
/// * `solver_settings` - Solver settings.
/// * `convergence_data` - Convergence data.
///
/// # Returns
///
/// * `true` if the value tolerance is satisfied, `false` otherwise.
///
/// # Note
///
/// If `solver_settings` has not defined `vtol`, then this function will just return `false`.
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
