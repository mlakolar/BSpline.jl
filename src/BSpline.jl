module BSpline

export
  bspline_basismatrix!, bspline_basismatrix,
  derivative_bspline_basismatrix!, derivative_bspline_basismatrix,
  derivative2_bspline_basismatrix!, derivative2_bspline_basismatrix,
  bspline_basis, derivative_bspline_basis, derivative2_bspline_basis,
  num_basis,
  BSplineT


type BSplineT{T<:FloatingPoint}
  M::Int64                  # B-spline order (2 for linear, 3 for quadratic, etc.)
  knots::Vector{T}          # knot vector
end

BSplineT{T<:FloatingPoint}(M::Int64, knots::Vector{T}) = BSplineT{T}(M, knots)
# BSplineT{T<:FloatingPoint}(M::Int64, knots::LinSpace{T}) = BSplineT{T}(M, [knots;])

function BSplineT{T<:FloatingPoint}(M::Int64, lb::T, ub::T, num_knots::Int64)
  h = (ub - lb) / (num_knots - 1)

  BSplineT{T}(M,
             vcat(lb .- h * [M-1:-1:1],
                  linspace(lb, ub, num_knots),
                  ub .+ h * [1:M-1]))
end


num_basis{T<:FloatingPoint}(bs::BSplineT{T}) = length(bs.knots) - bs.M

# B-spline basis function value B(j,n) at x.
#
# Input arguments:
# j:   interval index, 0 =< j < num_basis
#
# x:   value where the basis function is to be evaluated
#
# Output arguments:
# y:
#    B-spline basis function value, nonzero for a knot span of n

function bspline_basis{T<:FloatingPoint}(j::Integer, bs::BSplineT{T}, x::T)
  n = bs.M
  t = bs.knots
  y = zero(T)
  if n > 1
    dn = x - t[j+1]
    dd = t[j+n] - t[j+1]
    if dd != zero(T)  # indeterminate forms 0/0 are deemed to be zero
      b = bspline_basis(j, BSplineT(n-1,t), x)
      y = y + b*(dn/dd)
    end

    dn = t[j+n+1] - x
    dd = t[j+n+1] - t[j+2]
    if dd != zero(T)
      b = bspline_basis(j+1, BSplineT(n-1,t), x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] <= x < t[j+2] ? one(T) : zero(T)
  end
  y
end


# B-spline basis function value matrix B(n) for x.
#
# Input arguments:
# n:
#    B-spline order (2 for linear, 3 for quadratic, etc.)
# t:
#    knot vector
# x :
#    an m-dimensional vector of values where the basis function is to be
#    evaluated
#
# Output arguments:
# out:
#    a matrix of m rows and numel(t)-n columns

function bspline_basismatrix!{T<:FloatingPoint}(
    out::StridedMatrix{T}, bs::BSplineT{T}, x::Vector{T}
    )
  for j = 1 : num_basis(bs)
    for k = 1 : length(x)
      out[k, j] = bspline_basis(j-1, bs, x[k])
    end
  end
  out
end

bspline_basismatrix{T<:FloatingPoint}(bs::BSplineT{T}, x::T) =
  bspline_basismatrix!(Array(T, (1, num_basis(bs))), bs, [x;])

bspline_basismatrix{T<:FloatingPoint}(bs::BSplineT{T}, x::Vector{T}) =
  bspline_basismatrix!(Array(T, (length(x), num_basis(bs))), bs, x)

# derivatives

function derivative_bspline_basis{T<:FloatingPoint}(j::Integer, bs::BSplineT{T}, x::T)
  n = bs.M
  t = bs.knots

  y = zero(T)
  if n > 1
    dn = n - 1
    dd = t[j+n] - t[j+1]
    if dd != zero(T)  # indeterminate forms 0/0 are deemed to be zero
      b = bspline_basis(j, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end

    dn = 1 - n
    dd = t[j+n+1] - t[j+2]
    if dd != zero(T)
      b = bspline_basis(j+1, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] <= x < t[j+2] ? one(T) : zero(T)
  end
  y
end

function derivative_bspline_basismatrix!{T<:FloatingPoint}(
        out::StridedMatrix{T}, bs::BSplineT{T}, x::Vector{T}
        )
  for j = 1 : num_basis(bs)
    for k = 1 : length(x)
      out[k,j] = derivative_bspline_basis(j-1,bs,x[k])
    end
  end
  out
end

derivative_bspline_basismatrix{T<:FloatingPoint}(bs::BSplineT{T}, x::Vector{T}) =
  derivative_bspline_basismatrix!(Array(T, (length(x), num_basis(bs))), bs, x)

derivative_bspline_basismatrix{T<:FloatingPoint}(bs::BSplineT{T}, x::T) =
  derivative_bspline_basismatrix!(Array(T, (1, num_basis(bs))), bs, [x;])



# second derivative

function derivative2_bspline_basis{T<:FloatingPoint}(j::Integer, bs::BSplineT{T}, x::T)
  n = bs.M
  t = bs.knots
  y = zero(T)
  if n > 1
    dn = n - 1
    dd = t[j+n] - t[j+1]
    if dd != zero(T)  # indeterminate forms 0/0 are deemed to be zero
      b = derivative_bspline_basis(j, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end

    dn = 1 - n
    dd = t[j+n+1] - t[j+2]
    if dd != zero(T)
      b = derivative_bspline_basis(j+1, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] <= x < t[j+2] ? one(T) : zero(T)
  end
  y
end

function derivative2_bspline_basismatrix!{T<:FloatingPoint}(
            out::StridedMatrix{T}, bs::BSplineT{T}, x::Vector{T}
            )
  for j = 1 : num_basis(bs)
    for k = 1 : length(x)
      out[k,j] = derivative2_bspline_basis(j-1, bs, x[k])
    end
  end
end

derivative_bspline_basismatrix{T<:FloatingPoint}(bs::BSplineT{T}, x::Vector{T}) =
  derivative2_bspline_basismatrix!(Array(T, (length(x), num_basis(bs))), bs, x)

derivative_bspline_basismatrix{T<:FloatingPoint}(bs::BSplineT{T}, x::T) =
  derivative2_bspline_basismatrix!(Array(T, (1, num_basis(bs))), bs, [x;])

end

