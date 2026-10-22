"""
    BSpline

B-spline basis functions of any order and their first and second derivatives,
evaluated with the Cox–de Boor recursion.

Build a basis with [`BSplineT`](@ref), evaluate a single basis function with
[`bspline_basis`](@ref), or evaluate all of them at many points at once with
[`bspline_basismatrix`](@ref). The derivative versions add a `derivative_` or
`derivative2_` prefix.
"""
module BSpline

export
  BSplineT, num_basis,
  bspline_basis, derivative_bspline_basis, derivative2_bspline_basis,
  bspline_basismatrix!, bspline_basismatrix,
  derivative_bspline_basismatrix!, derivative_bspline_basismatrix,
  derivative2_bspline_basismatrix!, derivative2_bspline_basismatrix


"""
    BSplineT(M, knots)
    BSplineT(M, lb, ub, num_knots)

B-spline basis of order `M` on a nondecreasing knot vector `knots`. The order is
the degree plus one: `M = 1` gives piecewise constant functions, 2 linear,
3 quadratic and 4 cubic.

The basis has [`num_basis(bs)`](@ref num_basis)` = length(knots) - M`
functions. They are indexed from zero: basis function `j`, `0 ≤ j < num_basis`,
is a piecewise polynomial of degree `M - 1` that is positive on
`(knots[j+1], knots[j+M+1])` and zero outside `[knots[j+1], knots[j+M+1])`.

The second form spreads `num_knots` equally spaced knots over `[lb, ub]` and
adds `M - 1` more knots with the same spacing on each side. The resulting
`num_knots + M - 2` basis functions sum to one on `[lb, ub]` (on `[lb, ub)` when
`M = 1`).

The knots are copied and stored as floating point numbers: integer knots
become `Float64`, and `Float32` knots stay `Float32`.

# Examples
```julia
julia> bs = BSplineT(4, 0.0, 1.0, 5);   # cubic splines, knots at 0, 0.25, ..., 1

julia> num_basis(bs)
7

julia> bs.knots
11-element Vector{Float64}:
 -0.75
 -0.5
 -0.25
  0.0
  0.25
  0.5
  0.75
  1.0
  1.25
  1.5
  1.75
```
"""
struct BSplineT{T<:AbstractFloat}
  M::Int                    # B-spline order (2 for linear, 3 for quadratic, etc.)
  knots::Vector{T}          # knot vector, nondecreasing

  function BSplineT{T}(M::Integer, knots::AbstractVector{<:Real}) where {T<:AbstractFloat}
    M >= 1 || throw(ArgumentError("order M must be at least 1, got $M"))
    length(knots) > M ||
      throw(ArgumentError("order M = $M needs more than $M knots, got $(length(knots))"))
    issorted(knots) || throw(ArgumentError("knots must be nondecreasing"))
    new{T}(M, Vector{T}(knots))
  end
end

BSplineT(M::Integer, knots::AbstractVector{<:Real}) = BSplineT{float(eltype(knots))}(M, knots)

function BSplineT(M::Integer, lb::Real, ub::Real, num_knots::Integer)
  lb < ub || throw(ArgumentError("need lb < ub, got lb = $lb and ub = $ub"))
  num_knots >= 2 || throw(ArgumentError("need at least 2 knots, got $num_knots"))

  T = float(promote_type(typeof(lb), typeof(ub)))
  h = (T(ub) - T(lb)) / (num_knots - 1)
  knots = vcat(T(lb) .- h .* (M-1:-1:1),               # M - 1 knots below lb
               range(T(lb), T(ub); length=num_knots),
               T(ub) .+ h .* (1:M-1))                  # M - 1 knots above ub
  BSplineT{T}(M, knots)
end


"""
    num_basis(bs::BSplineT)

Number of basis functions in `bs`, `length(bs.knots) - bs.M`.
"""
num_basis(bs::BSplineT) = length(bs.knots) - bs.M


# `d`-th derivative at `x` of B(j,n), the order-`n` B-spline on the knots `t`
# with zero-based index `j`. B(j,n) vanishes outside [t[j+1], t[j+n+1]).
#
# Cox–de Boor recursion:
#   B(j,1)(x) = 1 if t[j+1] <= x < t[j+2], and 0 otherwise
#   B(j,n)(x) = (x - t[j+1]) / (t[j+n] - t[j+1]) * B(j,n-1)(x)
#             + (t[j+n+1] - x) / (t[j+n+1] - t[j+2]) * B(j+1,n-1)(x)
# Derivative, for d >= 1:
#   B(j,n)^(d)(x) = (n-1) * B(j,n-1)^(d-1)(x) / (t[j+n] - t[j+1])
#                 - (n-1) * B(j+1,n-1)^(d-1)(x) / (t[j+n+1] - t[j+2])
#
# A term whose denominator is zero (repeated knots) is dropped, which treats
# the indeterminate form 0/0 as zero. The recursion makes 2^(n-1) calls, which
# is cheap for the orders used in practice.
function _bspline(j::Integer, n::Integer, t::AbstractVector, x::Real, d::Integer)
  R = promote_type(eltype(t), typeof(x))

  # outside the support every term of the recursion is zero; checking first
  # skips the recursion for the many basis functions that vanish at x
  t[j+1] <= x < t[j+n+1] || return zero(R)
  if n == 1
    return d == 0 ? one(R) : zero(R)      # piecewise constant
  end

  dsub = d == 0 ? 0 : d - 1               # derivative order in the recursion
  y = zero(R)
  dd = t[j+n] - t[j+1]
  if !iszero(dd)
    dn = d == 0 ? x - t[j+1] : R(n - 1)
    y += dn / dd * _bspline(j, n-1, t, x, dsub)
  end
  dd = t[j+n+1] - t[j+2]
  if !iszero(dd)
    dn = d == 0 ? t[j+n+1] - x : R(1 - n)
    y += dn / dd * _bspline(j+1, n-1, t, x, dsub)
  end
  y
end

function _bspline(j::Integer, bs::BSplineT, x::Real, d::Integer)
  0 <= j < num_basis(bs) ||
    throw(ArgumentError("basis index j = $j is outside 0:$(num_basis(bs) - 1); indices are zero-based"))
  _bspline(j, bs.M, bs.knots, x, d)
end

# Fill out[k, j+1] with the `d`-th derivative of basis function j at x[k].
function _bspline_basismatrix!(out::AbstractMatrix, bs::BSplineT, x::AbstractVector{<:Real}, d::Integer)
  Base.require_one_based_indexing(out, x)
  size(out) == (length(x), num_basis(bs)) || throw(DimensionMismatch(
    "out has size $(size(out)), expected (length(x), num_basis(bs)) = $((length(x), num_basis(bs)))"))
  for j in axes(out, 2)
    for k in eachindex(x)
      out[k, j] = _bspline(j-1, bs.M, bs.knots, x[k], d)
    end
  end
  out
end

_basismatrix(bs::BSplineT{T}, x::AbstractVector{<:Real}) where {T} =
  Matrix{promote_type(T, eltype(x))}(undef, length(x), num_basis(bs))


"""
    bspline_basis(j, bs::BSplineT, x)

Value at `x` of basis function `j` of `bs`. The index is zero-based,
`0 ≤ j < num_basis(bs)`.

`x` can be any `Real`, including the dual numbers of automatic differentiation
packages. The result has type `promote_type(eltype(bs.knots), typeof(x))`.

Basis functions are right-continuous, and each one is zero at the right end of
its support. With a clamped knot vector, whose last knot is repeated `M` times,
every basis function is therefore zero at the last knot.

See also [`bspline_basismatrix`](@ref), [`derivative_bspline_basis`](@ref).

# Examples
```julia
julia> bs = BSplineT(4, 0.0, 1.0, 5);

julia> bspline_basis(2, bs, 0.25)   # cubic B-spline at its center
0.6666666666666666
```
"""
bspline_basis(j::Integer, bs::BSplineT, x::Real) = _bspline(j, bs, x, 0)

"""
    derivative_bspline_basis(j, bs::BSplineT, x)

First derivative at `x` of basis function `j` of `bs`, with `j` zero-based as in
[`bspline_basis`](@ref). At a knot where the derivative jumps, this is the
derivative from the right.
"""
derivative_bspline_basis(j::Integer, bs::BSplineT, x::Real) = _bspline(j, bs, x, 1)

"""
    derivative2_bspline_basis(j, bs::BSplineT, x)

Second derivative at `x` of basis function `j` of `bs`, with `j` zero-based as
in [`bspline_basis`](@ref). At a knot where the second derivative jumps, this is
the second derivative from the right. It is zero when `bs.M ≤ 2`.
"""
derivative2_bspline_basis(j::Integer, bs::BSplineT, x::Real) = _bspline(j, bs, x, 2)


"""
    bspline_basismatrix(bs::BSplineT, x::AbstractVector)
    bspline_basismatrix(bs::BSplineT, x::Real)

Matrix of all basis functions of `bs` evaluated at the points `x`: entry
`[k, j]` is [`bspline_basis(j-1, bs, x[k])`](@ref bspline_basis). The matrix
has `length(x)` rows and [`num_basis(bs)`](@ref num_basis) columns; a scalar `x`
gives a single row.

With `x` the observed values of a covariate, this is the design matrix of a
regression spline.

See also [`bspline_basismatrix!`](@ref), [`derivative_bspline_basismatrix`](@ref).

# Examples
```julia
julia> bs = BSplineT(2, 0.0, 1.0, 3);   # hat functions centered at 0, 0.5, 1

julia> bspline_basismatrix(bs, [0.0, 0.25, 1.0])
3×3 Matrix{Float64}:
 1.0  0.0  0.0
 0.5  0.5  0.0
 0.0  0.0  1.0
```
"""
bspline_basismatrix(bs::BSplineT, x::AbstractVector{<:Real}) =
  bspline_basismatrix!(_basismatrix(bs, x), bs, x)

bspline_basismatrix(bs::BSplineT, x::Real) = bspline_basismatrix(bs, [x])

"""
    bspline_basismatrix!(out, bs::BSplineT, x::AbstractVector)

In-place version of [`bspline_basismatrix`](@ref): fill `out`, of size
`length(x) × num_basis(bs)`, and return it.
"""
bspline_basismatrix!(out::AbstractMatrix, bs::BSplineT, x::AbstractVector{<:Real}) =
  _bspline_basismatrix!(out, bs, x, 0)


"""
    derivative_bspline_basismatrix(bs::BSplineT, x::AbstractVector)
    derivative_bspline_basismatrix(bs::BSplineT, x::Real)

Matrix of the first derivatives of all basis functions of `bs` at the points
`x`: entry `[k, j]` is
[`derivative_bspline_basis(j-1, bs, x[k])`](@ref derivative_bspline_basis).
The shape is as in [`bspline_basismatrix`](@ref).
"""
derivative_bspline_basismatrix(bs::BSplineT, x::AbstractVector{<:Real}) =
  derivative_bspline_basismatrix!(_basismatrix(bs, x), bs, x)

derivative_bspline_basismatrix(bs::BSplineT, x::Real) = derivative_bspline_basismatrix(bs, [x])

"""
    derivative_bspline_basismatrix!(out, bs::BSplineT, x::AbstractVector)

In-place version of [`derivative_bspline_basismatrix`](@ref): fill `out`, of
size `length(x) × num_basis(bs)`, and return it.
"""
derivative_bspline_basismatrix!(out::AbstractMatrix, bs::BSplineT, x::AbstractVector{<:Real}) =
  _bspline_basismatrix!(out, bs, x, 1)


"""
    derivative2_bspline_basismatrix(bs::BSplineT, x::AbstractVector)
    derivative2_bspline_basismatrix(bs::BSplineT, x::Real)

Matrix of the second derivatives of all basis functions of `bs` at the points
`x`: entry `[k, j]` is
[`derivative2_bspline_basis(j-1, bs, x[k])`](@ref derivative2_bspline_basis).
The shape is as in [`bspline_basismatrix`](@ref).

For a spline `f(x) = Σⱼ βⱼ Bⱼ(x)`, the roughness penalty `∫ f''(x)² dx` of
smoothing splines equals `β' Ω β` with `Ωᵢⱼ = ∫ Bᵢ''(x) Bⱼ''(x) dx`. If `D` is
this matrix on a fine grid with spacing `h`, then `Ω ≈ h * D' * D`.
"""
derivative2_bspline_basismatrix(bs::BSplineT, x::AbstractVector{<:Real}) =
  derivative2_bspline_basismatrix!(_basismatrix(bs, x), bs, x)

derivative2_bspline_basismatrix(bs::BSplineT, x::Real) = derivative2_bspline_basismatrix(bs, [x])

"""
    derivative2_bspline_basismatrix!(out, bs::BSplineT, x::AbstractVector)

In-place version of [`derivative2_bspline_basismatrix`](@ref): fill `out`, of
size `length(x) × num_basis(bs)`, and return it.
"""
derivative2_bspline_basismatrix!(out::AbstractMatrix, bs::BSplineT, x::AbstractVector{<:Real}) =
  _bspline_basismatrix!(out, bs, x, 2)

end
