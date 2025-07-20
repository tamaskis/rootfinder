# Changelog

## 0.4.0

1. Updated `rust` version to 2024.
1. Updated `once_cell` dependency from `1.20.3` to `1.21.3`.
1. Updated `numtest` dev dependency from `0.2.2` to `0.3.0`.

## 0.3.2

1. Updated `once_cell` dependency from `1.20.2` to `1.20.3`.

## 0.3.1

1. Updated `numtest` dev dependency from `0.2.1` to `0.2.2`.

## 0.3.0

1. Updated the `bisection_method` function to check for roots at the lower and upper bounds of the
initial interval.
1. Updated the `bisection_method` function to check whether the initial interval brackets a sign change.
1. Added relative bracket tolerance termination criteria.

## 0.2.1

1. Removed an unused, placeholder module.

## 0.2.0

1. Added Newton's method root solver.
1. Updated all functions to accept `f: &impl Fn(f64) -> f64` instead of `f: fn(f64) -> f64`.
1. Updated the `bisection_method` function to actually track the termination reason.
1. Updated `once_cell` dependency from `1.19.0` to `1.20.2`.
1. Updated `numtest` dev dependency from `0.1.6` to `0.2.1`.

## 0.1.0

1. Initial release with bisection method solvers and basic utilities.