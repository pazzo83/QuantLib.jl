type SampledCurve
  grid::Vector{Float64}
  values::Vector{Float64}
end

SampledCurve(gridSize::Int) = SampledCurve(zeros(gridSize), zeros(gridSize))

get_size(curve::SampledCurve) = length(curve.grid)

set_log_grid!(curve::SampledCurve, min::Float64, max::Float64) = set_grid!(curve, bounded_log_grid(min, max, get_size(curve) - 1))

set_grid!(curve::SampledCurve, grid::Vector{Float64}) = curve.grid = grid

function Math.sample!{T}(curve::SampledCurve, f::T)
  for i in eachindex(curve.grid)
    curve.values[i] = f(curve.grid[i])
  end

  return curve
end
