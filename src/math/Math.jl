# Math module
module Math

using Dierckx

# types for derivatives, etc
abstract FunctionType

type Derivative <: FunctionType end

type BernsteinPolynomial end

# misc function for comparison
function is_close{T <: Number}(x::T, y::T, n::Int = 42)
  if x == y
    return true
  end

  diff = abs(x - y)
  tol = n * eps(Float64) # machine epsilon

  if (x * y == 0.0) # x or y is 0
    return diff < (tol * tol)
  end

  return diff <= tol * abs(x) && diff <= tol * abs(y)
end

function close_enough{T <: Number}(x::T, y::T, n::Int = 42)
  if x == y
    return true
  end

  diff = abs(x -y)
  tol = n * eps()

  if x * y == 0
    return diff < (tol * tol)
  end

  return diff <= tol * abs(x) || diff <= tol * abs(y)
end

# misc functions - prob put in own file
function divide_array_by_self!{T, N <: Number}(a::Vector{T}, x::N)
  for i = 1:length(a)
    a[i] = a[i] / x
  end

  return a
end

function multiply_array_by_self!{T, N <: Number}(a::Vector{T}, x::N)
  for i = 1:length(a)
    a[i] = a[i] * x
  end

  return a
end

function get_factorial{I <: Integer}(i::I)
  if i > 20
    return Float64(factorial(BigInt(i)))
  else
    return factorial(i)
  end
end

function get_polynomial{I <: Integer}(::BernsteinPolynomial, i::I, n::I, x::Float64)
  coeff = get_factorial(n) / (get_factorial(n-1) * get_factorial(i))

  return coeff * (x ^ i) * (1.0 - x)^(n - i)
end

export is_close, close_enough, divide_array_by_self!, multiply_array_by_self!, get_factorial, get_polynomial, FunctionType, Derivative, BernsteinPolynomial

# Splines
type BSpline
  p::Integer
  n::Integer
  knots::Vector{Float64}
end

function spline_oper{I <: Integer}(spline::BSpline, i::I, x::Float64)
  i <= spline.n || error("i must not be greater than spline.n $i $(spline.n)")
  return N(spline, i, spline.p, x)
end

function N{I <: Integer}(spline::BSpline, i::I, p::I, x::Float64)
  if p == 0
    return (spline.knots[i] <= x && x < spline.knots[i + 1]) ? 1.0 : 0.0
  else
    return ((x - spline.knots[i]) / (spline.knots[i + p] - spline.knots[i])) * N(spline, i, p - 1, x) +
            ((spline.knots[i + p + 1] - x) / (spline.knots[i + p + 1] - spline.knots[i + 1])) * N(spline, i + 1, p - 1, x)
  end
end

export BSpline, spline_oper, N

# Constants
const EPS_VAL = eps()

# lmdif.jl
export lmdif!
include("lmdif.jl")

# lmdif2.jl
export lmdif2!
include("lmdif2.jl")

# distribution.jl
export distribution_derivative
include("distributions.jl")

# integral.jl
export Integrator, SegmentIntegral, operator, integrate
include("integral.jl")

# interpolation.jl
export Interpolation, LogInterpolation, BicubicSpline, NaturalCubicSpline, update!, locate, initialize!, value

include("interpolation.jl")

# solvers.jl
export Solver1D, BrentSolver, NewtonSolver, FiniteDifferenceNewtonSafe, solve

include("solvers.jl")

# optimization.jl
export Projection, CostFunction, Constraint, NoConstraint, PositiveConstraint, BoundaryConstraint, ProjectedConstraint, OptimizationMethod, LevenbergMarquardt, Simplex, Problem, EndCriteria,
project, minimize!, minimize_2!, include_params

include("optimization.jl")

end
