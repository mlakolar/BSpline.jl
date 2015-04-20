module BSpline

export
  derivative_bspline_basismatrix!,
  bspline_basismatrix!,
  derivative_bspline_basis,
  bspline_basis,
  num_basis,
  create_knots


type BSplineT
  M::Int64
  knots::Vector{Float64}
end


# B-spline basis function value B(j,n) at x.
#
# Input arguments:
# j:
#    interval index, 0 =< j < numel(t)-n
# n:
#    B-spline order (2 for linear, 3 for quadratic, etc.)
# t:
#    knot vector
# x :
#    value where the basis function is to be evaluated
#
# Output arguments:
# y:
#    B-spline basis function value, nonzero for a knot span of n


function bspline_basis(j::Integer, bs::BSplineT, x)

  n = bs.M
  t = bs.knots

  y = 0.

  if n > 1
    dn = x - t[j+1]
    dd = t[j+n] - t[j+1]
    if dd != 0  # indeterminate forms 0/0 are deemed to be zero
      b = bspline_basis(j, BSplineT(n-1,t), x)
      y = y + b*(dn/dd)
    end

    dn = t[j+n+1] - x
    dd = t[j+n+1] - t[j+2]
    if dd != 0
      b = bspline_basis(j+1, BSplineT(n-1,t), x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] <= x < t[j+2] ? 1. : 0.
  end

  return y
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

function bspline_basismatrix!(out::Matrix{Float64},
                              bs::BSplineT, x::Vector{Float64})
  for j = 1 : num_basis(bs)
    for k = 1 : length(x)
      out[k,j] = bspline_basis(j-1,bs,x[k])
    end
  end
end

function bspline_basismatrix!(out::Matrix{Float64},
                              bs::BSplineT, x::Float64)
  for j = 1 : num_basis(bs)
    out[1,j] = bspline_basis(j-1,bs,x)
  end
end


function num_basis(bs::BSplineT)
  length(bs.knots) - bs.M
end

function BSplineT(M::Int64, lb::Float64, ub::Float64, num_knots::Int64)
  h = (ub - lb) / (num_knots - 1)

  BSplineT(M,
           vcat(lb .- h * [M-1:-1:1],
                linspace(lb, ub, num_knots),
                ub .+ h * [1:M-1]))
end


function derivative_bspline_basis(j::Integer, bs::BSplineT, x)

  n = bs.M
  t = bs.knots

  y = 0.
  if n > 1
    dn = n - 1
    dd = t[j+n] - t[j+1]
    if dd != 0  # indeterminate forms 0/0 are deemed to be zero
      b = bspline_basis(j, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end

    dn = 1 - n
    dd = t[j+n+1] - t[j+2]
    if dd != 0
      b = bspline_basis(j+1, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] <= x < t[j+2] ? 1. : 0.
  end

  return y
end

function derivative_bspline_basismatrix!(out::Matrix{Float64},
                                         bs::BSplineT, x::Vector{Float64})
  for j = 1 : num_basis(bs)
    for k = 1 : length(x)
      out[k,j] = derivative_bspline_basis(j-1,bs,x[k]);
    end
  end
end

function derivative_bspline_basismatrix!(out::Matrix{Float64},
                                         bs::BSplineT, x::Float64)
  for j = 1 : num_basis(bs)
    out[1,j] = derivative_bspline_basis(j-1,bs,x);
  end
end

function derivative2_bspline_basis(j::Integer, bs::BSplineT, x)

  n = bs.M
  t = bs.knots

  y = 0.
  if n > 1
    dn = n - 1
    dd = t[j+n] - t[j+1]
    if dd != 0  # indeterminate forms 0/0 are deemed to be zero
      b = derivative_bspline_basis(j, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end

    dn = 1 - n
    dd = t[j+n+1] - t[j+2]
    if dd != 0
      b = derivative_bspline_basis(j+1, BSplineT(n-1, t), x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] <= x < t[j+2] ? 1. : 0.
  end

  return y
end

function derivative2_bspline_basismatrix!(out::Matrix{Float64},
                                         bs::BSplineT, x::Vector{Float64})
  for j = 1 : num_basis(bs)
    for k = 1 : length(x)
      out[k,j] = derivative2_bspline_basis(j-1,bs,x[k]);
    end
  end
end

function derivative2_bspline_basismatrix!(out::Matrix{Float64},
                                         bs::BSplineT, x::Float64)
  for j = 1 : num_basis(bs)
    out[1,j] = derivative2_bspline_basis(j-1,bs,x);
  end
end



end

