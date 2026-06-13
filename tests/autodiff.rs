use linalg_traits::Scalar;
use numdiff::{Dual, get_sderivative};
use numtest::*;
use rootfinder::{Interval, root_bisection, root_bisection_fast, root_newton, root_newton_fast};

// --------------
// Problem setup.
// --------------

// Define f(x,θ) = x² - θ.
fn f<S: Scalar>(x: S, theta: S) -> S {
    x.powi(2) - theta
}

// g(θ) = positive root of f(x,θ) = x² − θ with respect to x, for a fixed θ.
// dg(θ)/dθ = 1 / (2√(θ))

// g(θ) = positive root of f(x,θ) = x² − θ with respect to x, for a fixed θ.
// x² − θ = 0   -->    x* = √(θ)  -->   g(θ) = √(θ)
fn g_analyical<S: Scalar>(theta: S) -> S {
    theta.sqrt()
}

// g(θ) = positive root of f(x,θ) = x² − θ with respect to x, for a fixed θ.
// x² − θ = 0   -->    x* = √(θ)  -->   g(θ) = √(θ)
// dg(θ)/dθ = d/dθ(√(θ)) = 1 / (2√(θ))
fn dg_dtheta_analytical(theta: f64) -> f64 {
    1.0 / (2.0 * theta.sqrt())
}

// -----------------
// Bisection method.
// -----------------

fn g_root_bisection<S: Scalar>(theta: S, _p: &[f64]) -> S {
    // Define the root-finding problem h(x,θ) = 0 as a function of x alone.
    let h = |x: S| f(x, theta);

    // Return the root of h(x), which is g(θ).
    root_bisection(&h, Interval::new(0.0.into(), 5.0.into()), None, None).unwrap()
}

fn g_root_bisection_fast<S: Scalar>(theta: S, _p: &[f64]) -> S {
    // Define the root-finding problem h(x,θ) = 0 as a function of x alone.
    let h = |x: S| f(x, theta);

    // Return the root of h(x), which is g(θ).
    root_bisection_fast(&h, Interval::new(0.0.into(), 5.0.into()))
}

// ----------------
// Newton's method.
// ----------------

fn g_root_newton<S: Scalar>(theta: S, _p: &[f64]) -> S {
    // Define the root-finding problem h(x,θ) = 0 as a function of x alone.
    let h = |x: S| f(x, theta);

    // The derivative of h(x,θ) with respect to x is dh/dx = 2x.
    let dh = |x: S| x * 2.0;

    // Return the root of h(x), which is g(θ).
    root_newton(&h, &dh, 1.0.into(), None, None).unwrap()
}

fn g_root_newton_fast<S: Scalar>(theta: S, _p: &[f64]) -> S {
    // Define the root-finding problem h(x,θ) = 0 as a function of x alone.
    let h = |x: S| f(x, theta);

    // The derivative of h(x,θ) with respect to x is dh/dx = 2x.
    let dh = |x: S| x * 2.0;

    // Return the root of h(x), which is g(θ).
    root_newton_fast(&h, &dh, 1.0.into(), None)
}

// Get the derivatives.
get_sderivative!(g_root_bisection, dg_root_bisection);
get_sderivative!(g_root_bisection_fast, dg_root_bisection_fast);
get_sderivative!(g_root_newton, dg_root_newton);
get_sderivative!(g_root_newton_fast, dg_root_newton_fast);

#[test]
fn test_autodiff() {
    // Name of the root-finding function being tested.
    let function_names = [
        "root_bisection",
        "root_bisection_fast",
        "root_newton",
        "root_newton_fast",
    ];

    // Corresponding functions for g(θ) and dg(θ)/dθ, computed using different root-finding methods
    // and automatic differentiation.
    let g_functions = [
        g_root_bisection,
        g_root_bisection_fast,
        g_root_newton,
        g_root_newton_fast,
    ];
    let dg_functions = [
        dg_root_bisection,
        dg_root_bisection_fast,
        dg_root_newton,
        dg_root_newton_fast,
    ];

    // Iterate through the different root-finding methods and test that both g(θ) and dg(θ)/dθ match
    // their analytical expressions.
    for ((g, dg), name) in g_functions
        .iter()
        .zip(dg_functions.iter())
        .zip(function_names.iter())
    {
        // Test at multiple values of θ to make sure the test isn't just passing by coincidence at a
        // single value.
        for theta in [0.25, 1.0, 4.0] {
            // Status message to make it easier to debug if the test fails.
            println!(
                "Testing autodifferentiation of g(θ) with '{}' at θ = {}.",
                name, theta
            );

            // Test that g(θ) (obtained via a root-finding method) matches the analytical expression
            // for g(θ).
            assert_equal_to_decimal!(g(theta, &[]), g_analyical(theta), 16);

            // Test that dg(θ)/dθ (obtained by differentiating a root-finding method via
            //forward-mode automatic differentiation) matches the analytical expression for
            // dg(θ)/dθ.
            assert_equal_to_decimal!(dg(theta, &[]), dg_dtheta_analytical(theta), 16);
        }
    }
}
