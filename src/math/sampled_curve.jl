type SampledCurve
  grid::Vector{Float64}
  values::Vector{Float64}
end

SampledCurve(gridSize::Int) = SampledCurve(zeros(gridSize), zeros(gridSize))
