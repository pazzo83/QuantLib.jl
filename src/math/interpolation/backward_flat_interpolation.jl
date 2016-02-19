type BackwardFlatInterpolation <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  primitive::Vector{Float64}
end

function BackwardFlatInterpolation()
  x_vals = Vector{Float64}()
  y_vals = Vector{Float64}()
  primitive = Vector{Float64}()

  return BackwardFlatInterpolation(x_vals, y_vals, primitive)
end

# Backward Flat initialize #
function initialize!(interp::BackwardFlatInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})
  interp.x_vals = x_vals
  interp.y_vals = y_vals
  interp.primitive = zeros(length(y_vals))

  return interp
end

# BackwardFlatInterpolation update
function update!(interp::BackwardFlatInterpolation, idx::Int)
  interp.primitive[1] = 0.0
  for i = 2:idx
    dx = interp.x_vals[i] - interp.x_vals[i - 1]
    interp.primitive[i] = interp.primitive[i - 1] + dx * interp.y_vals[i]
  end

  return interp
end

function value(interp::BackwardFlatInterpolation, x::Float64)
  if x <= interp.x_vals[1]
    return interp.y_vals[1]
  end

  i = locate(interp, x)
  if x == interp.x_vals[i]
    return interp.y_vals[i]
  else
    return i >= length(interp.y_vals) ? 0.0 : interp.y_vals[i + 1]
  end
end

function get_primitive(interp::BackwardFlatInterpolation, x::Float64, ::Bool)
  i = locate(interp, x)
  # if i > 4
  #   println(x)
  #   println(interp.x_vals)
  #   println(interp.y_vals)
  # end
  dx = x - interp.x_vals[i]
  return interp.primitive[i] + dx * (i >= length(interp.y_vals) ? 0.0 : interp.y_vals[i + 1])
end
