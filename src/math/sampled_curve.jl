import Base.copy

type SampledCurve
  grid::Vector{Float64}
  values::Vector{Float64}
end

SampledCurve(gridSize::Int) = SampledCurve(zeros(gridSize), zeros(gridSize))
SampledCurve() = SampledCurve(Vector{Float64}(), Vector{Float64}())

get_size(curve::SampledCurve) = length(curve.grid)

set_log_grid!(curve::SampledCurve, min::Float64, max::Float64) = set_grid!(curve, bounded_log_grid(min, max, get_size(curve) - 1))

set_grid!(curve::SampledCurve, grid::Vector{Float64}) = curve.grid = grid

function Math.sample!{T}(curve::SampledCurve, f::T)
  # @simd for i in eachindex(curve.grid)
  #   @inbounds curve.values[i] = f(curve.grid[i])
  # end
  map!(f, curve.values, curve.grid)

  return curve
end

function value_at_center(curve::SampledCurve)
  n = get_size(curve)
  jmid = round(Int, floor(n / 2)) + 1
  if n % 2 == 1
    return curve.values[jmid]
  else
    return (curve.values[jmid] + curve.values[jmid-1]) / 2.0
  end
end

function first_derivative_at_center(curve::SampledCurve)
  n = get_size(curve)
  n >= 3 || error("the size of the curve must be at least 3")
  jmid = round(Int, floor(n / 2)) + 1
  if n % 2 == 1
    return (curve.values[jmid+1] - curve.values[jmid-1]) / (curve.grid[jmid+1] - curve.grid[jmid-1])
  else
    return (curve.values[jmid] - curve.values[jmid-1]) / (curve.grid[jmid] - curve.grid[jmid-1])
  end
end

function second_derivative_at_center(curve::SampledCurve)
  n = get_size(curve)
  n >= 4 || error("the size of the curve must be at least 4")
  jmid = round(Int, floor(n / 2)) + 1
  if n % 2 == 1
    deltaPlus = (curve.values[jmid+1] - curve.values[jmid]) / (curve.grid[jmid+1] - curve.grid[jmid])
    deltaMinus = (curve.values[jmid] - curve.values[jmid-1]) / (curve.grid[jmid] - curve.grid[jmid-1])
    dS = (curve.grid[jmid + 1] - curve.grid[jmid-1]) / 2.0
    return (deltaPlus - deltaMinus) / dS
  else
    deltaPlus = (curve.values[jmid+1] - curve.values[jmid-1]) / (curve.grid[jmid+1] - curve.grid[jmid-1])
    deltaMinus = (curve.values[jmid] - curve.values[jmid-2]) / (curve.grid[jmid] - curve.grid[jmid-2])
    return (deltaPlus - deltaMinus) / (curve.grid[jmid] - curve.grid[jmid-1])
  end
end

copy(curve::SampledCurve) = SampledCurve(copy(curve.grid), copy(curve.values))
