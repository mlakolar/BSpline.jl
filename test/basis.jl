module BasisTest

using Test
using BSpline
using ForwardDiff
using Random

# Greville abscissae: Σⱼ greville(bs, j) Bⱼ(x) = x wherever the basis sums to one
greville(bs, j) = sum(bs.knots[j+2:j+bs.M]) / (bs.M - 1)

@testset "constructors" begin
  bs = BSplineT(4, 0.0, 1.0, 5)
  @test bs.M == 4
  @test bs.knots ≈ -0.75:0.25:1.75
  @test num_basis(bs) == 7
  for M in 1:5, n in 2:6
    @test num_basis(BSplineT(M, -1.0, 2.0, n)) == n + M - 2
  end

  # element types
  @test BSplineT(4, 0, 1, 5).knots == bs.knots
  @test BSplineT(3, 0f0, 1f0, 4) isa BSplineT{Float32}
  @test BSplineT(3, 0:5) isa BSplineT{Float64}
  @test BSplineT(3, Float32[0, 1, 2, 3]) isa BSplineT{Float32}
  @test num_basis(BSplineT(2, 0:5)) == 4

  # the knots are copied
  t = collect(0.0:5.0)
  bs = BSplineT(3, t)
  t[1] = 10.0
  @test bs.knots[1] == 0.0

  @test_throws ArgumentError BSplineT(0, 0:5)
  @test_throws ArgumentError BSplineT(3, [0.0, 1.0, 2.0])
  @test_throws ArgumentError BSplineT(2, [0.0, 2.0, 1.0])
  @test_throws ArgumentError BSplineT(3, 1.0, 0.0, 5)
  @test_throws ArgumentError BSplineT(3, 0.0, 1.0, 1)
end

@testset "known values" begin
  # a single uniform cubic B-spline on [0, 4]
  bs = BSplineT(4, 0:4)
  @test num_basis(bs) == 1
  @test bspline_basis(0, bs, 0.5) ≈ 1/48
  @test bspline_basis(0, bs, 1.0) ≈ 1/6
  @test bspline_basis(0, bs, 2.0) ≈ 2/3
  @test bspline_basis(0, bs, 3.0) ≈ 1/6
  @test derivative_bspline_basis(0, bs, 1.0) ≈ 1/2
  @test derivative_bspline_basis(0, bs, 2.0) ≈ 0 atol=1e-14
  @test derivative2_bspline_basis(0, bs, 1.0) ≈ 1
  @test derivative2_bspline_basis(0, bs, 2.0) ≈ -2
  @test bspline_basis(0, bs, -0.1) == 0
  @test bspline_basis(0, bs, 4.0) == 0

  # a hat function on unequally spaced knots
  bs = BSplineT(2, [0, 1, 3])
  @test bspline_basis(0, bs, 0.5) ≈ 0.5
  @test bspline_basis(0, bs, 2.0) ≈ 0.5
  @test derivative_bspline_basis(0, bs, 0.5) ≈ 1
  @test derivative_bspline_basis(0, bs, 2.0) ≈ -0.5
end

@testset "partition of unity and reproduction of x" begin
  for M in 1:5
    bs = BSplineT(M, -1.0, 2.0, 7)
    x = range(-1.0, 2.0; length=301)
    x = M == 1 ? x[1:end-1] : x     # piecewise constant functions vanish at ub
    B = bspline_basismatrix(bs, x)
    @test all(>=(0), B)
    @test all(sum(B; dims=2) .≈ 1)
    if M >= 2
      @test B * [greville(bs, j) for j in 0:num_basis(bs)-1] ≈ x
    end

    # derivatives of a constant are zero; derivatives are from the right, and
    # to the right of ub the basis no longer sums to one, so leave out x = ub
    xl = range(-1.0, 2.0; length=301)[1:end-1]
    @test maximum(abs, sum(derivative_bspline_basismatrix(bs, xl); dims=2)) < 1e-10
    @test maximum(abs, sum(derivative2_bspline_basismatrix(bs, xl); dims=2)) < 1e-10
  end
end

@testset "repeated knots" begin
  # clamped cubic knots with a double interior knot at 2
  bs = BSplineT(4, [0, 0, 0, 0, 1, 2, 2, 3, 3, 3, 3])
  x = range(0.0, 3.0; length=301)[1:end-1]
  B = bspline_basismatrix(bs, x)
  @test all(isfinite, B)
  @test all(sum(B; dims=2) .≈ 1)
  @test B * [greville(bs, j) for j in 0:num_basis(bs)-1] ≈ x
  @test bspline_basis(0, bs, 0.0) == 1
  @test all(isfinite, derivative_bspline_basismatrix(bs, x))
  @test all(isfinite, derivative2_bspline_basismatrix(bs, x))

  # with a clamped knot vector every basis function vanishes at the last knot
  @test all(iszero, bspline_basismatrix(bs, 3.0))
end

@testset "derivatives match automatic differentiation" begin
  Random.seed!(1)
  bases = [[BSplineT(M, -1.0, 2.0, 6) for M in 1:5];
           BSplineT(4, [0, 0, 0, 0, 1, 2, 2, 3, 3, 3, 3]);
           BSplineT(3, [0.0, 0.1, 0.5, 0.6, 1.4, 2.0, 2.2])]
  for bs in bases
    nb = num_basis(bs)
    lo, hi = bs.knots[1], bs.knots[end]
    # random points and the knots themselves, where derivatives are from the right
    x = [lo .+ (hi - lo) .* rand(100); bs.knots[1:end-1]]

    D1 = derivative_bspline_basismatrix(bs, x)
    D2 = derivative2_bspline_basismatrix(bs, x)
    AD1 = [ForwardDiff.derivative(z -> bspline_basis(j, bs, z), xk) for xk in x, j in 0:nb-1]
    AD2 = [ForwardDiff.derivative(z -> ForwardDiff.derivative(w -> bspline_basis(j, bs, w), z), xk)
           for xk in x, j in 0:nb-1]
    @test D1 ≈ AD1 atol=1e-10
    @test D2 ≈ AD2 atol=1e-10
  end
end

@testset "derivatives of piecewise constant and linear functions" begin
  x = 0:0.1:4
  @test all(iszero, derivative_bspline_basismatrix(BSplineT(1, 0:4), x))
  @test all(iszero, derivative2_bspline_basismatrix(BSplineT(1, 0:4), x))
  @test all(iszero, derivative2_bspline_basismatrix(BSplineT(2, 0:4), x))
end

@testset "basis matrices" begin
  bs = BSplineT(4, 0.0, 1.0, 5)
  nb = num_basis(bs)
  xr = range(-1.5, 2.5; length=41)
  x = collect(xr)
  for (f, f!, g) in ((bspline_basismatrix, bspline_basismatrix!, bspline_basis),
                     (derivative_bspline_basismatrix, derivative_bspline_basismatrix!, derivative_bspline_basis),
                     (derivative2_bspline_basismatrix, derivative2_bspline_basismatrix!, derivative2_bspline_basis))
    B = f(bs, x)
    @test B isa Matrix{Float64}
    @test size(B) == (length(x), nb)
    @test B == [g(j, bs, xk) for xk in x, j in 0:nb-1]
    @test f(bs, xr) == B
    @test f(bs, x[20]) == B[20:20, :]

    out = fill(NaN, length(x), nb)
    @test f!(out, bs, x) === out
    @test out == B
    @test_throws DimensionMismatch f!(zeros(length(x), nb + 1), bs, x)
  end

  # outside the knots every basis function is zero
  @test all(iszero, bspline_basismatrix(bs, [-1.0, 2.0]))
end

@testset "argument types" begin
  bs32 = BSplineT(3, 0f0, 1f0, 4)
  @test bspline_basis(1, bs32, 0.5f0) isa Float32
  @test bspline_basismatrix(bs32, Float32[0.1, 0.2]) isa Matrix{Float32}
  @test bspline_basismatrix(bs32, [0.1, 0.2]) isa Matrix{Float64}

  bs = BSplineT(3, 0, 4, 5)
  @test bspline_basis(1, bs, 1) isa Float64
  @test bspline_basis(1, bs, 1) == bspline_basis(1, bs, 1.0)
  @test bspline_basismatrix(bs, 0:4) isa Matrix{Float64}

  @test_throws ArgumentError bspline_basis(-1, bs, 0.5)
  @test_throws ArgumentError derivative_bspline_basis(num_basis(bs), bs, 0.5)
end

end
