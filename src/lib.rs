//! [![github]](https://github.com/tamaskis/rootfinder)&ensp;[![crates-io]](https://crates.io/crates/rootfinder)&ensp;[![docs-rs]](https://docs.rs/rootfinder)
//!
//! [github]: https://img.shields.io/badge/github-8da0cb?style=for-the-badge&labelColor=555555&logo=github
//! [crates-io]: https://img.shields.io/badge/crates.io-fc8d62?style=for-the-badge&labelColor=555555&logo=rust
//! [docs-rs]: https://img.shields.io/badge/docs.rs-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs
//!
//! Root-finding methods for both univariate, scalar-valued functions and multivariate, vector-valued functions.
//!
//! # Root-finding functions for univariate, scalar-valued functions
//!
//! `rootfinder` provides the following functions for finding the root of a univariate,
//! scalar-valued function, $f:\mathbb{R}\to\mathbb{R}$:
//!
//! | Method | Functions |
//! | ------ | --------- |
//! | bisection | [`root_bisection`]<BR>[`root_bisection_fast`] |
//! | Newton's | [`root_newton`]<BR>[`root_newton_fast`] |

// Linter setup.
#![warn(missing_docs)]

// Module declarations.
pub(crate) mod univariate_roots;
pub(crate) mod utils;

// Re-exports.
pub use crate::univariate_roots::bisection::{root_bisection, root_bisection_fast};
pub use crate::univariate_roots::newton::{
    root_newton, root_newton_fast, DEFAULT_NEWTON_SOLVER_SETTINGS,
};
pub use crate::utils::bracketing::Interval;
pub use crate::utils::enums::SolverError;
pub use crate::utils::solver_settings::{SolverSettings, DEFAULT_SOLVER_SETTINGS};
