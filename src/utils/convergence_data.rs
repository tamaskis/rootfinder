use crate::utils::enums::TerminationReason;

#[derive(Debug, Default)]
pub struct ConvergenceData {
    /// Solutions at all iterations.
    pub x_all: Vec<f64>,

    /// Bracketing interval lower bounds at all iterations.
    pub a_all: Vec<f64>,

    /// Bracketing interval upper bounds at all iterations.
    pub b_all: Vec<f64>,

    /// Function evaluations at all iterations.
    pub f_all: Vec<f64>,

    /// Derivative evaluations at all iterations.
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
