# BSpline.jl

> **Archived.** This package is no longer maintained. Use the registered
> [BSplineKit.jl](https://github.com/jipolanco/BSplineKit.jl), which computes the
> same basis functions and derivatives. For `x` in `[lb, ub]`, the basis
> `BSplineT(M, lb, ub, num_knots)` and its matrices translate as follows, with
> one-based basis indices:
>
> ```julia
> using BSplineKit
> h = (ub - lb) / (num_knots - 1)
> knots = vcat(lb .- h .* (M-1:-1:1), range(lb, ub; length=num_knots), ub .+ h .* (1:M-1))
> B = BSplineBasis(BSplineOrder(M), knots; augment=Val(false))
> X = collocation_matrix(B, x, Derivative(0), Matrix{Float64})  # bspline_basismatrix(bs, x)
> D = collocation_matrix(B, x, Derivative(2), Matrix{Float64})  # derivative2_bspline_basismatrix(bs, x)
> Ω = galerkin_matrix(B, (Derivative(2), Derivative(2)))        # exact roughness penalty on [lb, ub]
> ```
>
> For a knot vector `t`, use `BSplineBasis(BSplineOrder(M), t; augment=Val(false))`.
> BSplineKit takes the left limit at the last knot, so a clamped basis also sums
> to one there.

[![CI](https://github.com/mlakolar/BSpline.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/mlakolar/BSpline.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![codecov](https://codecov.io/gh/mlakolar/BSpline.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mlakolar/BSpline.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

B-spline basis functions of any order and their first and second derivatives,
evaluated with the Cox–de Boor recursion. Use it to build the design matrix of
a regression spline and the roughness penalty of a smoothing spline.

## Installation

Requires Julia 1.10 or later. The package is not registered; install it from
GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/mlakolar/BSpline.jl")
```

## Example: fitting a regression spline

```julia
using BSpline

n = 200
x = sort(rand(n))
y = sin.(2π .* x) .+ 0.3 .* randn(n)

bs = BSplineT(4, 0.0, 1.0, 20)    # cubic splines, 20 equally spaced knots on [0, 1]
B = bspline_basismatrix(bs, x)    # 200 × 22 design matrix
β = B \ y                         # least squares fit

xg = range(0, 1; length=101)
f̂ = bspline_basismatrix(bs, xg) * β              # fitted curve
f̂′ = derivative_bspline_basismatrix(bs, xg) * β  # its derivative
```

With 22 basis functions for 200 points this fit is wiggly. A smoothing spline
penalizes the roughness `∫ f''(x)² dx = β'Ωβ`, where
`Ωᵢⱼ = ∫ Bᵢ''(x) Bⱼ''(x) dx`. A Riemann sum on a fine grid gives `Ω` from the
second-derivative matrix:

```julia
grid = range(0, 1; length=1001)[1:end-1]
D = derivative2_bspline_basismatrix(bs, grid)
Ω = step(grid) * D' * D

λ = 1e-4
β = (B' * B + λ * Ω) \ (B' * y)
```

## The basis

`BSplineT(M, knots)` is the basis of order `M` on a nondecreasing knot vector
`t = knots`. The order is the degree plus one: `M = 1` gives piecewise constant
functions, 2 linear, 3 quadratic and 4 cubic. There are
`num_basis(bs) = length(t) - M` basis functions, defined by

```
B(j,1)(x) = 1 if t[j+1] ≤ x < t[j+2], and 0 otherwise
B(j,M)(x) = (x - t[j+1]) / (t[j+M] - t[j+1]) * B(j,M-1)(x)
          + (t[j+M+1] - x) / (t[j+M+1] - t[j+2]) * B(j+1,M-1)(x)
```

where a term with a zero denominator (from repeated knots) is dropped. Basis
function `j` is a piecewise polynomial of degree `M - 1` that is zero outside
`[t[j+1], t[j+M+1])`, so at any `x` at most `M` basis functions are nonzero.

`BSplineT(M, lb, ub, num_knots)` places `num_knots` equally spaced knots on
`[lb, ub]` and `M - 1` more on each side with the same spacing. The resulting
`num_knots + M - 2` basis functions sum to one on `[lb, ub]`. This is usually
the basis you want for a covariate observed on `[lb, ub]`.

## Functions

| Function | Returns |
|---|---|
| `BSplineT(M, knots)`, `BSplineT(M, lb, ub, num_knots)` | basis of order `M` |
| `num_basis(bs)` | number of basis functions |
| `bspline_basis(j, bs, x)` | value of basis function `j` at `x` |
| `derivative_bspline_basis(j, bs, x)` | its first derivative |
| `derivative2_bspline_basis(j, bs, x)` | its second derivative |
| `bspline_basismatrix(bs, x)` | `length(x) × num_basis(bs)` matrix of all basis functions at the points `x` |
| `derivative_bspline_basismatrix(bs, x)` | the same for first derivatives |
| `derivative2_bspline_basismatrix(bs, x)` | the same for second derivatives |

Each matrix function has an in-place version ending in `!`, for example
`bspline_basismatrix!(out, bs, x)`. Use `?name` in the Julia REPL for the full
docstring.

## Conventions

- **Basis indices are zero-based**: `bspline_basis(j, bs, x)` takes
  `0 ≤ j < num_basis(bs)`, and column `j + 1` of `bspline_basismatrix` holds
  basis function `j`. This matches the 2015 version of the package.
- **Right-continuity**: each basis function is zero at the right end of its
  support, and at a knot the derivatives are taken from the right. With a
  clamped knot vector, whose last knot is repeated `M` times, every basis
  function is zero at the last knot; with `BSplineT(M, lb, ub, num_knots)` and
  `M ≥ 2` the basis still sums to one at `ub`.
- **Number types**: knots are stored as floating point numbers, and `x` can be
  any `Real`. In particular, the scalar functions accept the dual numbers of
  [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl), so they work
  inside automatically differentiated code.

## Example notebook

[Using BSplines.ipynb](Using%20BSplines.ipynb) plots a cubic basis with its
first and second derivatives and checks the derivatives against ForwardDiff.

## License

MIT, see [LICENSE](LICENSE).
