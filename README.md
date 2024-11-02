# rootfinder

[<img alt="github" src="https://img.shields.io/badge/github-tamaskis/rootfinder-8da0cb?style=for-the-badge&labelColor=555555&logo=github" height="20">](https://github.com/tamaskis/rootfinder)
[<img alt="crates.io" src="https://img.shields.io/crates/v/rootfinder.svg?style=for-the-badge&color=fc8d62&logo=rust" height="20">](https://crates.io/crates/rootfinder)
[<img alt="docs.rs" src="https://img.shields.io/badge/docs.rs-rootfinder-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs" height="20">](https://docs.rs/rootfinder)

Root-finding methods for both univariate, scalar-valued functions and multivariate, vector-valued functions.

## Documentation

Please see https://docs.rs/rootfinder.

## Examples

### Finding a root using the bisection method.

```rust
use rootfinder::{root_bisection, Interval};

// Define the function f(x) = x² - 1.
let f = |x: f64| x.powi(2) - 1.0;

// We want the root in the interval [0,∞). Therefore, we use an initial interval of
// [a,b] = [0,9999999]. Finding this root using the bisection method,
let result = root_bisection(f, Interval::new(0.0, 9999999.0), None, None);
let root = result.unwrap();

// `root` is `0.9999999999999999`, which is very close to the true root of 1.
```

#### License

<sup>
Licensed under either of <a href="LICENSE-APACHE">Apache License, Version 2.0</a> or 
<a href="LICENSE-MIT">MIT license</a> at your option.
</sup>

<br>

<sub>
Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in
this crate by you, as defined in the Apache-2.0 license, shall be dual licensed as above, without
any additional terms or conditions.
</sub>