# Math module
module Math

using Dierckx

# types for derivatives, etc
abstract type FunctionType end

struct Derivative <: FunctionType end

struct BernsteinPolynomial end

# misc function for comparison
function is_close(x::T, y::T, n::Int = 42) where {T <: Number}
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

function close_enough(x::T, y::T, n::Int = 42) where {T <: Number}
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
function divide_array_by_self!(a::Vector{T}, x::N) where {T, N <: Number}
  for i = 1:length(a)
    a[i] = a[i] / x
  end

  return a
end

function multiply_array_by_self!(a::Vector{T}, x::N) where {T, N <: Number}
  for i = 1:length(a)
    a[i] = a[i] * x
  end

  return a
end

function get_factorial(i::Int)
  if i > 20
    return Float64(factorial(BigInt(i)))
  else
    return factorial(i)
  end
end

function get_polynomial(::BernsteinPolynomial, i::Int, n::Int, x::Float64)
  coeff = get_factorial(n) / (get_factorial(n-1) * get_factorial(i))

  return coeff * (x ^ i) * (1.0 - x)^(n - i)
end

export is_close, close_enough, divide_array_by_self!, multiply_array_by_self!, get_factorial, get_polynomial, FunctionType, Derivative, BernsteinPolynomial

# Splines
mutable struct BSpline
  p::Integer
  n::Integer
  knots::Vector{Float64}
end

function spline_oper(spline::BSpline, i::Int, x::Float64)
  i <= spline.n || error("i must not be greater than spline.n $i $(spline.n)")
  return N(spline, i, spline.p, x)
end

function N(spline::BSpline, i::Int, p::Int, x::Float64)
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

# utilities.jl
include("utilities.jl")

# grid.jl
export bounded_log_grid
include("grid.jl")

# lmdif2.jl
export lmdif2!
include("lmdif2.jl")

# distribution.jl
export distribution_derivative, peizer_pratt_method_2_inversion
include("distributions.jl")

# integral.jl
export Integrator, IntegrationFunction, SegmentIntegral, operator, integrate, GaussLaguerreIntegration, get_order
include("integral.jl")

# svd.jl
# export SVD
# include("svd.jl")

# tridiagonal_operator.jl
export TridiagonalOperator, TridiagIdentity, set_mid_row!, set_last_row!, set_first_row!, apply_to, solve_for!, solve_for
include("tridiagonal_operator.jl")

# matrix.jl
export NoneSalvagingAlgo, rank_reduced_sqrt, OrthogonalProjection, get_vector
include("matrix.jl")

# transformed_grid.jl
export TransformedGrid, LogGrid
include("transformed_grid.jl")

# rng.jl
export AbstractRandomSequenceGenerator, PseudoRandomRSG, InverseCumulativeRSG, SobolRSG, SobolInverseCumulativeRSG, next_sequence!, last_sequence, init_sequence_generator!
include("rng.jl")

# statistics.jl
export NonWeightedStatistics, GenericRiskStatistics, RiskStatistics, gen_RiskStatistics, GenericSequenceStats, reset!, add_sample!, adding_data!, error_estimate, stats_mean, stats_std_deviation,
stats_skewness, stats_kurtosis, stats_covariance
include("statistics.jl")

# sampled_curve.jl
export SampledCurve, get_size, set_log_grid!, set_grid!, sample!, value_at_center, first_derivative_at_center, second_derivative_at_center
include("sampled_curve.jl")

# general_linear_least_squares.jl
export GeneralLinearLeastSquares, get_coefficients
include("general_linear_least_squares.jl")

# interpolation.jl
export Interpolation, Spline, Lagrange, LogInterpolation, LogLinear, BicubicSpline, NaturalCubicSpline, CubicInterpolation, BackwardFlatInterpolation, update!, locate, initialize!, value, get_primitive,
derivative

include("interpolation/interpolation.jl")
include("interpolation/linear_interpolation.jl")
include("interpolation/log_interpolation.jl")
include("interpolation/backward_flat_interpolation.jl")
include("interpolation/spline_interpolation.jl")

# solvers.jl
export Solver1D, BrentSolver, NewtonSolver, FiniteDifferenceNewtonSafe, solve
include("solvers/solver.jl")
include("solvers/brent_solver.jl")
include("solvers/finite_difference.jl")
include("solvers/newton_solver.jl")

# optimization.jl
export lmdif!, Projection, CostFunction, Constraint, NoConstraint, PositiveConstraint, BoundaryConstraint, ProjectedConstraint, OptimizationMethod, LevenbergMarquardt, Simplex, Problem, EndCriteria,
project, minimize!, minimize_2!, include_params
include("optimization/lmdif.jl")
include("optimization/optimization.jl")
include("optimization/problem.jl")
include("optimization/simplex.jl")
include("optimization/levenberg_marquardt.jl")

end
