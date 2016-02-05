abstract Integrator

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
