abstract Integrator
abstract GaussianQuadrature
abstract GaussianOrthogonalPolynomial

type GaussLaguerrePolynomial <: GaussianOrthogonalPolynomial
  s::Float64
end

get_alpha(poly::GaussLaguerrePolynomial, i::Int) = 2 * i + 1 + poly.s
get_beta(poly::GaussLaguerrePolynomial, i::Int) = i * (i + poly.s)
get_mu_0(poly::GaussLaguerrePolynomial) = exp(lgamma(poly.s + 1.0))
get_w(poly::GaussLaguerrePolynomial, x::Float64) = ^(x, poly.s) * exp(-x)

type SegmentIntegral{I <: Integer} <: Integrator
  absoluteAccuracy::Float64
  absoluteError::Float64
  maxEvals::I
  evals::I
  intervals::I
end

SegmentIntegral{I <: Integer}(intervals::I) = SegmentIntegral{I}(1.0, 0.0, 1, 0, intervals)


function operator{I <: Integrator}(integrator::I, f::Function, a::Float64, b::Float64)
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

function Math.integrate(integrator::SegmentIntegral, f::Function, a::Float64, b::Float64)
  dx = (b - a) / integrator.intervals
  sum_ = 0.5 * (f(a) + f(b))
  end_ = b - 0.5 * dx
  x = a + dx
  while x < end_
    sum_ += f(x)
    x += dx
  end

  return sum_ * dx
end

type GaussLaguerreIntegration <: GaussianQuadrature
  x::Vector{Float64}
  w::Vector{Float64}
end

function GaussLaguerreIntegration(n::Int, s::Float64 = 0.0)
  x, w = build_gaussian_quadrature(n, GaussLaguerrePolynomial(s))

  return GaussLaguerreIntegration(x, w)
end

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
