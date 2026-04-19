# EML Algebraic Coordinates

## Goal

`eml(x, y) = exp(x) - ln(y)` is expressive enough to generate the full backend, but it is not pleasant to read or write directly. This experiment adds a small surface language with a coordinate-like constructor:

```text
a + ε b
```

The intent is analogous to writing `x + i y` for complex numbers: the user gets a compact two-part notation, while the compiler lowers it into the existing backend.

## Core Idea

`a + ε b` is a frontend constructor for:

```text
eml(a, b)
```

More precisely, the lowering path is:

```text
EmlCoordExpr -> ScientificExpr -> EmlExpr
```

The internal `EmlExpr` representation is unchanged.

## Surface AST

The implementation introduces two new frontend-only ASTs:

- `SurfaceScalarExpr`: a small parsed scalar language that lowers to `ScientificExpr`
- `EmlCoordExpr`: the coordinate surface language

`EmlCoordExpr` supports:

- scalar lifting
- coordinate construction with `a + ε b` and `a - ε b`
- addition and subtraction
- scalar multiplication
- experimental general multiplication

## Parsing

The parser is intentionally minimal and self-contained. It supports:

- variables such as `x`, `velocity`, `u_1`
- scalar constants `1`, `e`, `i`, `pi`
- scalar functions already present in the backend, such as `exp`, `ln`, `sqrt`, `sin`, `cos`, `tan`, and inverse/hyperbolic variants
- `ε`, `eps`, or `epsilon` as the coordinate marker

This is not a full mathematical parser for the crate as a whole. It exists only for this experimental surface feature.

## Lowering Rules

- `Scalar(s)` lowers to `s.to_scientific()`
- `Coord(a, b)` lowers to `exp(a) - ln(b)`
- `Add`, `Sub`, and `ScalarMul` lower through the matching `ScientificExpr` operations
- `Mul` lowers to ordinary `ScientificExpr::mul`

## Multiplication Status

Multiplication is only partially supported.

The parser and AST accept:

```text
(x + ε y) * (u + ε v)
```

and lower it correctly into an ordinary product of the lowered scientific expressions. However, the result is not normalized back into a single coordinate literal. In other words, multiplication is supported as a surface-language operation, but not yet as a closed algebra over the coordinate form itself.

## Non-Goals

- no change to `EmlExpr`
- no change to the numerical evaluator semantics
- no promise that coordinate expressions form a canonical or closed algebraic system
- no attempt to replace the rest of the library API with parsed surface expressions
