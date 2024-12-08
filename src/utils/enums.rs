/// Solver errors.
#[derive(Debug)]
pub enum SolverError {
    /// An error that occurs when a bracketing interval is not found.
    BracketingIntervalNotFound,
}

impl std::fmt::Display for SolverError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SolverError::BracketingIntervalNotFound => {
                write!(
                    f,
                    "No interval was found that brackets a sign change of the function."
                )
            }
        }
    }
}

/// Solver termination reasons.
///
/// # Note
///
/// We derive the `PartialEq` trait to facilitate easier unit testing.
#[derive(Debug, Default, PartialEq)]
pub enum TerminationReason {
    /// Solver not yet terminated.
    #[default]
    NotYetTerminated,

    /// Solver terminated on reaching the maximum number of iterations.
    MaxIterationsReached,

    /// Solver terminated on reaching the maximum number of function evaluations.
    MaxFunctionEvaluationsReached,

    /// Solver terminated on satisfying the absolute bracket tolerance.
    AbsoluteBracketToleranceSatisfied,

    /// Solver terminated on satisfying the absolute step tolerance.
    AbsoluteStepToleranceSatisfied,

    /// Solver terminated on satisfying the value tolerance.
    ValueToleranceSatisfied,

    /// Solver terminated due to a zero derivative.
    ///
    /// # Note
    ///
    /// This is only used for Newton's method.
    ZeroDerivative,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solver_error_bracketing_interval_not_found() {
        assert_eq!(
            format!("{}", SolverError::BracketingIntervalNotFound),
            "No interval was found that brackets a sign change of the function."
        );
    }

    #[test]
    fn test_termination_reason_default() {
        let termination_reason = TerminationReason::default();
        assert!(matches!(
            termination_reason,
            TerminationReason::NotYetTerminated
        ));
    }
}
