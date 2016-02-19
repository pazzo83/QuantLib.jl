type LinearInterpolation <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  s::Vector{Float64}
end

# Linear initialize
function initialize!(interp::LinearInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})
  interp.x_vals = x_vals
  interp.y_vals = y_vals
  interp.s = zeros(length(y_vals))

  return interp
end

# Linear Interpolation update
function update!{I <: Integer}(interp::LinearInterpolation, idx::I)
  for i = 2:idx
    @inbounds dx = interp.x_vals[i] - interp.x_vals[i - 1]
    @inbounds interp.s[i - 1] = (interp.y_vals[i] - interp.y_vals[i - 1]) / dx
  end

  return interp
end

function value(interp::LinearInterpolation, val::Float64)
  i = locate(interp, val)
  # println("I is ", i)
  # println("Val is ", val)
  # println("Y vals ", interp.y_vals)
  # println("X vals ", interp.x_vals)
  # println("S vals ", interp.s)
  return interp.y_vals[i] + (val - interp.x_vals[i]) * interp.s[i]
end
