# emlmath

`emlmath` is a small Rust crate that explores the idea from *All elementary functions from a single binary operator* by Andrzej Odrzywołek. The project implements the paper's `eml(x, y) = exp(x) - log(y)` construction and builds a library of elementary functions on top of that single binary operator.

The repository currently contains:

- A `Complex` number type with the operations needed to evaluate elementary functions.
- A `ComplexBall` type for midpoint-plus-radius error bounds during tree evaluation.
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

For bound tracking, the crate also exposes ball-arithmetic evaluators:

- `eval_ball` / `eval_ball_scientific` for explicit `ComplexBall` assignments
- `eval_real_ball` / `eval_real_ball_scientific` for real inputs with scalar radii

Supported derived operations include:

- Arithmetic: `add`, `sub`, `mul`, `div`, `pow`, `reciprocal`
- Constants: `0`, `1`, `2`, `1/2`, `e`, `i`, `pi`
- Elementary functions: `exp`, `ln`, `sqrt`
- Trigonometric and hyperbolic functions: `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`
- Inverse functions: `asin`, `acos`, `atan`, `asinh`, `acosh`, `atanh`

## Ball Arithmetic

`emlmath` can evaluate expression trees with conservative midpoint-plus-radius bounds using `ComplexBall`. This is useful when you want an enclosure for the result instead of a single point value.

Example:

```rust
use emlmath::{ComplexBall, ScientificExpr, eval_ball_scientific, sx};

fn main() {
    let expr = ScientificExpr::ln(ScientificExpr::add(sx(), ScientificExpr::one()));
    let value = eval_ball_scientific(
        &expr,
        &[("x", ComplexBall::from_real(0.5, 1e-6))],
    )
    .unwrap();

    println!("value = {}", value);
}
```

The ball evaluator is conservative and explicit about unsafe regions:

- `ln` and the `eml` operator reject balls that contain zero.
- `ln` and `eml_log` reject balls that cross the nonpositive real axis.
- Division rejects denominator balls that contain zero.

When those cases occur, evaluation returns a `BallEvalError` instead of silently producing an invalid bound.

## Branch Handling

EML compilation is branch-sensitive because `eml(x, y)` is defined in terms of `exp` and `ln`, and the compiled `ln` tree crosses the negative real axis internally.

This crate implements the paper's partial branch-handling solution in two places:

- `EmlExpr` evaluation uses a mirrored logarithm branch internally so that compiled `ln(z)` matches the standard principal branch across the negative real axis.
- The compiler uses a sign-corrected internal `i` when building derived `EmlExpr` formulas that depend on the branch convention.

These choices are enough to make compiled `ln` and the mirrored-log branch behavior consistent with the intended principal-branch output.

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
The bundled demo binary builds the expression as `ScientificExpr`, compiles it to `EmlExpr` for inspection, and prints both the direct point evaluation and a ball-arithmetic enclosure:

```bash
cargo run
```

At the time of writing, that example produces an `EmlExpr` with `715` nodes for `sin(x) + sqrt(y)`, together with the direct value and a conservative error ball. This shows both how quickly the single-operator representation expands and how the new bound-tracking path can be used alongside the direct evaluator.

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
- The paper's branch workaround is only partial: compiled `ln` is validated, but some higher derived `EmlExpr` constructors still do not match direct `ScientificExpr` evaluation when complex intermediate values appear.
- In particular, compiled trigonometric/hyperbolic and some inverse-function formulas still have known equivalence gaps; the test suite keeps those checks as ignored regression targets until the compiler identities are re-derived correctly.
- The generated `EmlExpr` trees can become very large very quickly.
- Ball-arithmetic bounds are conservative and can widen substantially on deep compiled `EmlExpr` trees.
- Some expressions that are well-defined at a point may still be rejected in ball mode if the enclosing ball touches a singularity or branch cut.

## License

This project is available under the terms in [`LICENSE`](/home/tim/Documents/GitHub/emlmath/LICENSE).
