# using FloatFloat

type FittingMethodCommons{T, I <: Integer}
  solution::Vector{T}
  guessSolution::Vector{T}
  numberOfIterations::I
  minimumCostValue::Float64
  weights::Vector{T}
  costFunction::FittingCost
end

function FittingMethodCommons{I <: Integer}(size::I, gsize::I)
  solution = zeros(size)
  # solution = Vector{DD}(size)
  guessSolution = zeros(gsize)
  # guessSolution = Vector{DD}(gsize)
  numberOfIterations = 0
  minimumCostValue = 0.0
  weights = zeros(size)
  # weights = Vector{DD}(size)
  curve = NullCurve()
  costFunction = FittingCost(size, curve)

  return FittingMethodCommons(solution, guessSolution, numberOfIterations, minimumCostValue, weights, costFunction)
end

type ExponentialSplinesFitting{I <: Integer} <: FittingMethod
  constrainAtZero::Bool
  size::I
  commons::FittingMethodCommons
end

function ExponentialSplinesFitting{I <: Integer}(constrainAtZero::Bool, size::I)
  if constrainAtZero
    gsize = 9
  else
    gsize = 10
  end

  commons = FittingMethodCommons(size, gsize)

  return ExponentialSplinesFitting(constrainAtZero, gsize, commons)
end

type SimplePolynomialFitting{I <: Integer} <: FittingMethod
  constrainAtZero::Bool
  degree::I
  size::I
  commons::FittingMethodCommons
end

function SimplePolynomialFitting{I <: Integer}(constrainAtZero::Bool, degree::I, size::I)
  if constrainAtZero
    gsize = degree
  else
    gsize = degree + 1
  end

  commons = FittingMethodCommons(size, gsize)

  return SimplePolynomialFitting(constrainAtZero, degree, gsize, commons)
end

type NelsonSiegelFitting{I <: Integer} <: FittingMethod
  constrainAtZero::Bool
  size::I
  commons::FittingMethodCommons
end

function NelsonSiegelFitting{I <: Integer}(size::I)
  constrainAtZero = true
  gsize = 4

  commons = FittingMethodCommons(size, gsize)

  return NelsonSiegelFitting(constrainAtZero, gsize, commons)
end

type SvenssonFitting{I <: Integer} <: FittingMethod
  constrainAtZero::Bool
  size::I
  commons::FittingMethodCommons
end

function SvenssonFitting{I <: Integer}(size::I)
  constrainAtZero = true
  gsize = 6

  commons = FittingMethodCommons(size, gsize)

  return SvenssonFitting(constrainAtZero, gsize, commons)
end

type CubicBSplinesFitting{I <: Integer} <: FittingMethod
  constrainAtZero::Bool
  size::I
  knots::Vector{Float64}
  splines::BSpline
  N::I
  commons::FittingMethodCommons

  function CubicBSplinesFitting(constrainAtZero::Bool, knots::Vector{Float64}, size::I)
    m = length(knots)
    m >= 8 || error("At least 8 knots are required")

    splines = BSpline(3, m - 4, knots)
    basis_functions = m - 4
    if constrainAtZero
      gsize = basis_functions - 1

      N = 2
      abs(spline_oper(splines, N, 0.0) > JQuantLib.Math.EPS_VAL) || error("N_th cubic B-spline must be nonzero at t=0")
    else
      gsize = basis_functions
      N = 1
    end

    commons = FittingMethodCommons(size, gsize)

    new(constrainAtZero, gsize, knots, splines, N, commons)
  end
end

CubicBSplinesFitting{I <: Integer}(constrainAtZero::Bool, knots::Vector{Float64}, size::I) = CubicBSplinesFitting{I}(constrainAtZero, knots, size)


# function ExponentialSplinesFitting(constrainAtZero::Bool, size::Integer)
#   solution = zeros(size)
#   if constrainAtZero
#     gsize = 9
#   else
#     gsize = 10
#   end
#   guessSolution = zeros(gsize)
#   numberOfIterations = 0
#   minimumCostValue = 0.0
#   weights = zeros(size)
#   curve = NullCurve()
#   costFunction = FittingCost(size, curve)
#
#   return ExponentialSplinesFitting(constrainAtZero, gsize,
#           FittingMethodCommons(solution, guessSolution, numberOfIterations, minimumCostValue, weights, costFunction))
# end

guess_size(fitting::ExponentialSplinesFitting) = fitting.constrainAtZero ? 9 : 10
guess_size(fitting::SimplePolynomialFitting) = fitting.constrainAtZero ? fitting.degree : fitting.degree + 1

# Discount functions
function discount_function{T}(method::ExponentialSplinesFitting, x::Vector{T}, t::Float64)
  d = 0.0
  N = guess_size(method)
  kappa = x[N]
  coeff = 0.0
  if !method.constrainAtZero
    for i = 1:N
      @inbounds d += x[i] * exp(-kappa * (i) * t)
    end
  else
    for i = 1:N - 1
      @inbounds d += x[i]  * exp(-kappa * (i + 1) * t)
      @inbounds coeff += x[i]
    end
    coeff = 1.0 - coeff
    d += coeff * exp(-kappa * t)
  end
  return d
end

function discount_function{T}(method::SimplePolynomialFitting, x::Vector{T}, t::Float64)
  d = 0.0
  N = method.size

  if !method.constrainAtZero
    for i = 1:N
      @inbounds d += x[i] * get_polynomial(BernsteinPolynomial(), i-1, i-1, t)
    end
  else
    d = 1.0
    for i = 1:N
      @inbounds d += x[i] * get_polynomial(BernsteinPolynomial(), i, i, t)
    end
  end

  return d
end

function discount_function{T}(method::NelsonSiegelFitting, x::Vector{T}, t::Float64)
  kappa = x[method.size]
  @inbounds zero_rate = x[1] + (x[2] + x[3]) * (1.0 - exp(-kappa * t)) / ((kappa + JQuantLib.Math.EPS_VAL) * (t + JQuantLib.Math.EPS_VAL)) - (x[3]) * exp(-kappa * t)
  d = exp(-zero_rate * t)

  return d
end

function discount_function{T}(method::SvenssonFitting, x::Vector{T}, t::Float64)
  kappa = x[method.size - 1]
  kappa_1 = x[method.size]
  eps_v = JQuantLib.Math.EPS_VAL

  zero_rate = x[1] + (x[2] + x[3]) * (1.0 - exp(-kappa * t)) / ((kappa + eps_v) * (t + eps_v)) - (x[3]) * exp(-kappa * t) + x[4] * (((1.0 - exp(-kappa_1 * t)) / ((kappa_1 + eps_v) * (t + eps_v))) - exp(-kappa_1 * t))

  d = exp(-zero_rate * t)
  return d
end

function discount_function{T}(method::CubicBSplinesFitting, x::Vector{T}, t::Float64)
  d = 0.0
  if !method.constrainAtZero
    for i = 1:method.size
      @inbounds d += x[i] * spline_oper(method.splines, i, t)
    end
  else
    t_star = 0.0
    sum = 0.0
    for i = 1:method.size
      if i < method.N
        @inbounds d += x[i] * spline_oper(method.splines, i, t)
        @inbounds sum += x[i] * spline_oper(method.splines, i, t_star)
      else
        @inbounds d += x[i] * spline_oper(method.splines, i + 1, t)
        @inbounds sum += x[i] * spline_oper(method.splines, i + 1, t_star)
      end
    end
    coeff = 1.0 - sum
    coeff /= spline_oper(method.splines, method.N, t_star)
    d += coeff * spline_oper(method.splines, method.N, t)
  end

  return d
end
