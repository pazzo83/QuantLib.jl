function bounded_log_grid(xMin::Float64, xMax::Float64, steps::Int)
  result = zeros(steps + 1)
  gridLogSpacing = (log(xMax) - log(xMin)) / steps
  edx = exp(gridLogSpacing)
  result[1] = xMin
  for j = 2:steps+1
    result[j] = result[j-1] * edx
  end

  return result
end
