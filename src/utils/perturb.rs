use linalg_traits::Scalar;

/// Perturb a real number.
///
/// # Arguments
///
/// * `x` - Real number, $x\in\mathbb{R}$.
///
/// # Returns
///
/// `x` slightly perturbed.
///
/// # Note
///
/// $x$ is pertubed as $x+100\varepsilon(1+\|x\|)$ where $\varepsilon$ is defined by
/// [`f64::EPSILON`].
pub fn perturb_real<S: Scalar>(x: S) -> S {
    x + S::new(100.0) * S::new(f64::EPSILON) * (S::one() + x.abs())
}

#[cfg(test)]
mod tests {
    use super::*;
    use numtest::*;

    #[test]
    fn test_perturb_real_positive() {
        assert_equal!(perturb_real(2.0), 2.0000000000000666);
    }

    #[test]
    fn test_perturb_real_negative() {
        assert_equal!(perturb_real(-2.0), -1.9999999999999334);
    }

    #[test]
    fn test_perturb_real_zero() {
        assert_equal!(perturb_real(0.0), 2.220446049250313e-14);
    }
}
