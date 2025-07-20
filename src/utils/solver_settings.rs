use core::f64;
use once_cell::sync::Lazy;

/// Default solver settings.
pub static DEFAULT_SOLVER_SETTINGS: Lazy<SolverSettings> = Lazy::new(SolverSettings::default);

/// Solver settings.
#[derive(Default, Clone)]
pub struct SolverSettings {
    /// Absolute step tolerance.
    ///
    /// The absolute step tolerance is satisfied when
    ///
    /// `(x_next - x_current).abs() <= xatol`
    ///
    /// where `x_current` and `x_next` are the current and next iterates, respectively.
    pub xatol: Option<f64>,

    /// Value tolerance.
    ///
    /// The value tolerance is satisfied when
    ///
    /// `v_current.abs() <= vtol`
    ///
    /// where `v_current` is the function evaluation at the current iterate.
    pub vtol: Option<f64>,

    /// Absolute bracket tolerance.
    ///
    /// The absolute bracket tolerance is satisfied when
    ///
    /// `(b - a).abs() <= batol`
    ///
    /// where `a` and `b` are the lower and upper bounds, respectively, of the current bracketing
    /// interval.
    pub batol: Option<f64>,

    /// Relative bracket tolerance.
    ///
    /// The relative bracket tolerance is satisfied when
    ///
    /// `(b - a).abs() <= brtol * (a.abs().max(b.abs()))`
    ///
    /// where `a` and `b` are the lower and upper bounds, respectively, of the current bracketing
    /// interval.
    pub brtol: Option<f64>,

    /// Maximum number of solver iterations allowed.
    pub max_iter: Option<u32>,

    /// Maximum number of function evaluations allowed.
    pub max_feval: Option<u32>,

    /// Maximum number of derivative evaluations allowed.
    pub max_deval: Option<u32>,

    /// `true` if the initial interval should be updated to ensure a sign change, `false` otherwise.
    pub rebracket: Option<bool>,

    /// Maximum number of bracketing iterations allowed.
    pub max_bracket_iter: Option<u32>,
}

impl SolverSettings {
    /// Constructor.
    ///
    /// # Arguments
    /// * `xatol` - Absolute step tolerance. The absolute step tolerance is satisfied when
    ///   `(x_next - x_current).abs() <= xatol`, where `x_current` and `x_next` are the current and
    ///   next iterates, respectively.
    /// * `vtol` - Value tolerance. The value tolerance is satisfied when `v_current.abs() <= vtol`,
    ///   where `v_current` is the function evaluation at the current iterate.
    /// * `batol` - Absolute bracket tolerance. The absolute bracket tolerance is satisfied when
    ///   `(b - a).abs() <= batol`, where `a` and `b` are the lower and upper bounds, respectively,
    ///   of the current bracketing interval.
    /// * `brtol` - Relative bracket tolerance. The relative bracket tolerance is satisfied when
    ///   `(b - a).abs() <= brtol * (a.abs().max(b.abs()))`, where `a` and `b` are the lower and
    ///   upper bounds, respectively, of the current bracketing interval.
    /// * `max_iter` - Maximum number of solver iterations allowed.
    /// * `max_feval` - Maximum number of function evaluations allowed.
    /// * `max_deval` - Maximum number of derivative evaluations allowed.
    /// * `rebracket` - `true` if the initial interval should be updated to ensure a sign change,
    ///   `false` otherwise.
    /// * `max_bracket_iter` - Maximum number of iterations to find a bracketing interval allowed.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        xatol: Option<f64>,
        vtol: Option<f64>,
        batol: Option<f64>,
        brtol: Option<f64>,
        max_iter: Option<u32>,
        max_feval: Option<u32>,
        max_deval: Option<u32>,
        rebracket: Option<bool>,
        max_bracket_iter: Option<u32>,
    ) -> Self {
        Self {
            xatol,
            vtol,
            batol,
            brtol,
            max_iter,
            max_feval,
            max_deval,
            rebracket,
            max_bracket_iter,
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::utils::solver_settings::DEFAULT_SOLVER_SETTINGS;

    #[test]
    fn test_default_solver_settings() {
        assert!(DEFAULT_SOLVER_SETTINGS.xatol.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.vtol.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.batol.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.max_iter.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.max_feval.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.max_deval.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.rebracket.is_none());
        assert!(DEFAULT_SOLVER_SETTINGS.max_bracket_iter.is_none());
    }
}
