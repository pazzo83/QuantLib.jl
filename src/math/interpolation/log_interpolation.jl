type LogInterpolation{IN <: Interpolation} <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  interpolator::IN
end

typealias LogLinearInterpolation LogInterpolation{LinearInterpolation}

function LogLinear(x_vals::Vector{Float64}, y_vals::Vector{Float64})
  # build log of y values, defaulting to 0 for initial state
  n = length(x_vals)
  log_y_vals = zeros(n)
  s = zeros(n) # first derivative

  # initialize the linear interpolator
  interpolator = LinearInterpolation(x_vals, log_y_vals, s)

  return LogInterpolation(x_vals, y_vals, interpolator)
end

# if no values are provided
function LogLinear()
  x_vals = Vector{Float64}()
  y_vals = Vector{Float64}()
  s = Vector{Float64}()

  interpolator = LinearInterpolation(x_vals, y_vals, s)

  return LogInterpolation(x_vals, y_vals, interpolator)
end

# Log initialize
function initialize!(interp::LogInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})
  interp.x_vals = x_vals
  interp.y_vals = y_vals

  log_y = zeros(length(y_vals))
  for i = 1:length(y_vals)
    @inbounds log_y[i] = log(y_vals[i])
  end

  initialize!(interp.interpolator, x_vals, log_y)

  return interp
end

# Log Interpolation update
function update!{I <: Integer}(interp::LogInterpolation, idx::I)
  # first get the log of the y values
  for i = 1:idx
    @inbounds interp.interpolator.y_vals[i] = log(interp.y_vals[i])
  end
  # use these log values to update the linear interpolator
  update!(interp.interpolator, idx)

  return interp
end

value(interp::LogInterpolation, val::Float64) = exp(value(interp.interpolator, val))
