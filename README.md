# emlmath

`emlmath` is a small Rust crate that explores the idea from *All elementary functions from a single binary operator* by Andrzej Odrzywołek. The project implements the paper's `eml(x, y) = exp(x) - log(y)` construction and builds a library of elementary functions on top of that single binary operator.

The repository currently contains:

- A `Complex` number type with the operations needed to evaluate elementary functions.
- A low-level `EmlExpr` tree that only uses `1`, variables, and the `eml(left, right)` operator.
- A higher-level `ScientificExpr` tree with familiar functions like `sin`, `sqrt`, `atanh`, and `pow`, plus a compiler to `EmlExpr`.
- A small demo binary in [`src/main.rs`](/home/tim/Documents/GitHub/emlmath/src/main.rs).
- A copy of the reference paper in [`doc/2603.21852v2.pdf`](/home/tim/Documents/GitHub/emlmath/doc/2603.21852v2.pdf).

## Why this exists

Most math libraries start with many primitive operations and functions. This project goes in the opposite direction: it treats `eml` as the only primitive operator and derives everything else from it. That makes the crate useful as:

- A reference implementation of the paper's formulas.
- A playground for inspecting how large elementary expressions become when reduced to a single operator.
- A testbed for branch-cut and complex-analysis behavior when translating standard formulas into `eml` form.

## Current API

The crate exposes two main expression types:

- `ScientificExpr`: ergonomic constructors for standard math expressions.
- `EmlExpr`: the compiled single-operator form.

Typical workflow:

1. Build a `ScientificExpr` or `EmlExpr`.
2. Optionally convert a `ScientificExpr` into `EmlExpr` with `to_eml()`.
3. Evaluate the expression with real variable assignments using `eval_real` / `eval_real_scientific`, or evaluate directly with complex assignments via `eval`.

Supported derived operations include:

- Arithmetic: `add`, `sub`, `mul`, `div`, `pow`, `reciprocal`
- Constants: `0`, `1`, `2`, `1/2`, `e`, `i`, `pi`
- Elementary functions: `exp`, `ln`, `sqrt`
- Trigonometric and hyperbolic functions: `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`
- Inverse functions: `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`

## Example

```rust
use emlmath::{ScientificExpr, eval_real_scientific, sx, sy};

fn main() {
    let expr = ScientificExpr::add(ScientificExpr::sin(sx()), ScientificExpr::sqrt(sy()));
    let eml = expr.to_eml();
    let value = eval_real_scientific(&expr, &[("x", 0.5), ("y", 9.0)]).unwrap();

    println!("scientific = {:?}", expr);
    println!("eml nodes   = {}", eml.node_count());
    println!("value       = {}", value);
}
```

The bundled demo binary uses the same idea with `EmlExpr` directly:

```bash
cargo run
```

At the time of writing, that example produces an `EmlExpr` with `715` nodes for `sin(x) + sqrt(y)`, which shows how quickly the single-operator representation expands.

## Development

Build and test with standard Cargo commands:

```bash
cargo run
cargo test
```

The crate currently has no external dependencies.

## Project status

This is an experimental math/code exploration project, not a production numerical library. A few practical constraints are worth knowing:

- The public API is still minimal and may change freely.
- Numerical behavior follows the crate's custom `Complex` implementation rather than a battle-tested external library.
- Branch handling is intentionally tailored to make the compiled `eml` forms agree with the project's expected logarithm behavior.
- The generated `EmlExpr` trees can become very large very quickly.

## License

This project is available under the terms in [`LICENSE`](/home/tim/Documents/GitHub/emlmath/LICENSE).
