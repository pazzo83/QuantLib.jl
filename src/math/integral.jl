using SpecialFunctions

abstract type Integrator end
abstract type GaussianQuadrature end
abstract type GaussianOrthogonalPolynomial end
abstract type IntegrationFunction end

struct GaussLaguerrePolynomial <: GaussianOrthogonalPolynomial
  s::Float64
end

get_alpha(poly::GaussLaguerrePolynomial, i::Int) = 2 * i + 1 + poly.s
get_beta(poly::GaussLaguerrePolynomial, i::Int) = i * (i + poly.s)
get_mu_0(poly::GaussLaguerrePolynomial) = exp(lgamma(poly.s + 1.0))
get_w(poly::GaussLaguerrePolynomial, x::Float64) = ^(x, poly.s) * exp(-x)

mutable struct SegmentIntegral <: Integrator
  absoluteAccuracy::Float64
  absoluteError::Float64
  maxEvals::Int
  evals::Int
  intervals::Int
end

SegmentIntegral(intervals::Int) = SegmentIntegral(1.0, 0.0, 1, 0, intervals)


function operator(integrator::Integrator, f::Function, a::Float64, b::Float64)
  integrator.evals = 0
  if a == b
    return 0.0
  end

  if (b > a)
    return integrate(integrator, f, a, b)
  else
    return -integrate(integrator, f, b, a)
  end
end

function (integrator::SegmentIntegral)(f::Function, a::Float64, b::Float64)
  integrator.evals = 0
  if a == b
    return 0.0
  end

  if (b > a)
    return integrate(integrator, f, a, b)
  else
    return -integrate(integrator, f, b, a)
  end
end

function Math.integrate(integrator::SegmentIntegral, f::Union{Function, IntegrationFunction}, a::Float64, b::Float64)
  dx = (b - a) / (integrator.intervals * 1.0)
  sum_ = 0.5 * (f(a) + f(b))
  end_ = b - 0.5 * dx
  x = a + dx
  while x < end_
    sum_ += f(x)
    x += dx
  end

  return sum_ * dx
end

struct GaussLaguerreIntegration <: GaussianQuadrature
  x::Vector{Float64}
  w::Vector{Float64}
end

function (gli::GaussLaguerreIntegration)(f::IntegrationFunction)
  sum_ = 0.0
  @inbounds @simd for i = length(gli.x):-1:1
    sum_ += gli.w[i] * f(gli.x[i])
  end

  return sum_
end

function GaussLaguerreIntegration(n::Int, s::Float64 = 0.0)
  x, w = build_gaussian_quadrature(n, GaussLaguerrePolynomial(s))

  return GaussLaguerreIntegration(x, w)
end

get_order(integration::GaussianQuadrature) = length(integration.x)

function build_gaussian_quadrature(n::Int, poly::GaussianOrthogonalPolynomial)
  x = zeros(n)
  w = zeros(n)

  e = zeros(n - 1)

  @simd for i = 2:n
    @inbounds x[i] = get_alpha(poly, i - 1)
    @inbounds e[i - 1] = sqrt(get_beta(poly, i - 1))
  end

  x[1] = get_alpha(poly, 0)

  # sym tridiagonal
  # tridiag = SymTridiagonal(x, e)
  # eigtri = eigfact(tridiag)

  tqr = TqrEigenDecomposition(x, e, OnlyFirstRowEigenVector(), Overrelaxation())
  x_new = tqr.d
  ev = tqr.ev

  mu_0 = get_mu_0(poly)

  @simd for i = 1:n
    @inbounds w[i] = mu_0 * ev[i, 1] * ev[i, 1] / get_w(poly, x_new[i])
  end

  return x_new, w
end
