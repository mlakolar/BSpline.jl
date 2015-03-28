module BSpline

export
  derivative_bspline_basismatrix!,
  bspline_basismatrix!,
  derivative_bspline_basis,
  bspline_basis,
  num_basis,
  create_knots



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


function bspline_basis(j::Integer, n::Integer, t::Vector{Float64}, x::Float64)

  y = 0.

  if n > 1
    dn = x - t[j+1]
    dd = t[j+n] - t[j+1]
    if dd != 0  # indeterminate forms 0/0 are deemed to be zero
      b = bspline_basis(j,n-1,t,x)
      y = y + b*(dn/dd)
    end

    dn = t[j+n+1] - x
    dd = t[j+n+1] - t[j+2]
    if dd != 0
      b = bspline_basis(j+1,n-1,t,x)
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
                              n::Integer, t::Vector{Float64}, x::Vector{Float64})
  for j = 0 : length(t)-n-1
    for k = 1 : length(x)
      out[k,j+1] = bspline_basis(j,n,t,x[k])
    end
  end
end

function bspline_basismatrix!(out::Matrix{Float64},
                              n::Integer, t::Vector{Float64}, x::Float64)
  for j = 0 : length(t)-n-1
    out[1,j+1] = bspline_basis(j,n,t,x)
  end
end


function num_basis(n::Integer, t::Vector{Float64})
  length(t) - n
end

function create_knots(n::Integer, range::(Real,Real), num_knots::Integer)
  h = (range[2] - range[1]) / (num_knots - 1)
  vcat(range[1] .- h * [n-1:-1:1],
       linspace(range[1],range[2],num_knots),
       range[2] .+ h * [1:n-1])
end


function derivative_bspline_basis(j::Integer, n::Integer, t::Vector{Float64}, x::Float64)

  y = 0.
  if n > 1
    dn = n - 1
    dd = t[j+n] - t[j+1]
    if dd != 0  # indeterminate forms 0/0 are deemed to be zero
      b = bspline_basis(j,n-1,t,x)
      y = y + b*(dn/dd)
    end

    dn = 1 - n
    dd = t[j+n+1] - t[j+2]
    if dd != 0
      b = bspline_basis(j+1,n-1,t,x)
      y = y + b*(dn/dd)
    end
  else
    y = t[j+1] .<= x < t[j+2] ? 1. : 0.
  end

  return y
end


function derivative_bspline_basismatrix!(out::Matrix{Float64},
                                        n::Integer, t::Vector{Float64}, x::Vector{Float64})
  for j = 0 : length(t)-n-1
    for k = 1 : length(x)
      out[k,j+1] = derivative_bspline_basis(j,n,t,x[k]);
    end
  end
end

function derivative_bspline_basismatrix!(out::Matrix{Float64},
                                        n::Integer, t::Vector{Float64}, x::Float64)
  for j = 0 : length(t)-n-1
    out[1,j+1] = derivative_bspline_basis(j,n,t,x);
  end
end


end

