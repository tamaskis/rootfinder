use crate::utils::enums::TerminationReason;

/// Convergence data.
///
/// # Note
///
/// * The 0th element of the storage vectors (e.g. `x_all`, `f_all`, etc.) store values
///   corresponding to the initial guess (i.e. the 0th iteration).
/// * The number of iterations is always one less than the number of iterates. This is because the
///   1st iteration takes us from iterate 0 (i.e. the initial guess) to iterate 1, the 2nd iteration
///   takes us from iterate 1 to iterate 2, etc.
#[derive(Debug, Default)]
pub struct ConvergenceData {
    /// Solutions at all iterations.
    pub x_all: Vec<f64>,

    /// Bracketing interval lower bounds at all iterations.
    pub a_all: Vec<f64>,

    /// Bracketing interval upper bounds at all iterations.
    pub b_all: Vec<f64>,

    /// Function evaluations at all iterations.
    ///
    /// # Warning
    ///
    /// If the function was not evaluated at the initial guess, `f_all[0]` will be [`f64::NAN`].
    pub f_all: Vec<f64>,

    /// Derivative evaluations at all iterations.
    ///
    /// # Warning
    ///
    /// If the derivative was not evaluated at the initial guess, `df_all[0]` will be [`f64::NAN`].
    pub df_all: Vec<f64>,

    /// Number of iterations to find a bracketing interval.
    pub n_bracket_iter: u32,

    /// Number of solver iterations.
    pub n_iter: u32,

    /// Number of function evaluations.
    pub n_feval: u32,

    /// Number of derivative evaluations.
    pub n_deval: u32,

    /// Solver termination reason.
    pub termination_reason: TerminationReason,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default() {
        let convergence_data = ConvergenceData::default();
        assert_eq!(convergence_data.x_all, Vec::new());
        assert_eq!(convergence_data.a_all, Vec::new());
        assert_eq!(convergence_data.b_all, Vec::new());
        assert_eq!(convergence_data.f_all, Vec::new());
        assert_eq!(convergence_data.df_all, Vec::new());
        assert_eq!(convergence_data.n_bracket_iter, 0);
        assert_eq!(convergence_data.n_iter, 0);
        assert_eq!(convergence_data.n_feval, 0);
        assert_eq!(convergence_data.n_deval, 0);
    }
}
