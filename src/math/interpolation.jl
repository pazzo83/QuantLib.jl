# interpolation
using Dierckx

abstract Interpolation
abstract Interpolation2D <: Interpolation

type LinearInterpolation <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  s::Vector{Float64}
end

# Log Linear interpolation
type LogInterpolation <: Interpolation
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  interpolator::LinearInterpolation
end

type BicubicSpline
  spline::Dierckx.Spline2D
end

BicubicSpline{T <: Real}(x::Vector{T}, y::Vector{T}, z::Matrix{T}) = BicubicSpline(Dierckx.Spline2D(x, y, z))

type NaturalCubicSpline{T <: Number} <: Interpolation
  x_vert::Vector{T}
  y_vert::Vector{T}
  b::Vector{T}
  c::Vector{T}
  d::Vector{T}
end

# adapted from here: http://sepwww.stanford.edu/sep/sergey/128A/answers6.pdf
function NaturalCubicSpline{T <: Real}(x_vert::Vector{T}, y_vert::Vector{T})
  n = length(x_vert)
  h = zeros(n - 1)
  b = zeros(n - 1)
  a = zeros(n - 2)
  g = zeros(n - 2)
  c = zeros(n)
  d = zeros(n - 1)

  for i = 1:n - 1
    h[i] = x_vert[i + 1] - x_vert[i]
    b[i] = (y_vert[i + 1] - y_vert[i]) / h[i]
  end

  for i = 1:n - 2
    a[i] = 2.0 * (h[i] + h[i + 1])
    # u[i] = 6 * (b[i + 1] - b[i])
    g[i] = b[i + 1] - b[i]
  end

  Alu = lufact(Tridiagonal(h[2:end-1], a[1:end], h[2:end-1]))
  c[2:end - 1] = Alu \ g

  for i=1:n - 1
    d[i] = (c[i+1] - c[i])/h[i]
    b[i] -= (2.0 * c[i] + c[i + 1]) * h[i]
    c[i] *= 3.0
  end

  return NaturalCubicSpline(x_vert, y_vert, b, c, d)
end

function call{T <: Real}(spl::NaturalCubicSpline, my_x::T)
  # get coefficients
  # b, c, d = gen_splines(x_vert, y_vert)

  x_idx = searchsortedlast(spl.x_vert, my_x)

  if x_idx == 0
    x_idx = 1
  elseif x_idx == length(spl.x_vert)
    x_idx = length(spl.x_vert) - 1
  end


  diff = my_x - spl.x_vert[x_idx]

  return spl.y_vert[x_idx] + diff*(spl.b[x_idx] + diff * (spl.c[x_idx] + (diff * spl.d[x_idx])))
end

function LogInterpolation(x_vals::Vector{Float64}, y_vals::Vector{Float64})
  # build log of y values, defaulting to 0 for initial state
  n = length(x_vals)
  log_y_vals = zeros(n)
  s = zeros(n) # first derivative

  # initialize the linear interpolator
  interpolator = LinearInterpolation(x_vals, log_y_vals, s)

  return LogInterpolation(x_vals, y_vals, interpolator)
end

# if no values are provided
function LogInterpolation()
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

# Linear initialize
function initialize!(interp::LinearInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})
  interp.x_vals = x_vals
  interp.y_vals = y_vals
  interp.s = zeros(length(y_vals))

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

# update if value passed in
function update!{I <: Integer}(interp::LogInterpolation, idx::I, val::Float64)
  interp.y_vals[idx] = val

  update!(interp, idx)

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

# locate x
function locate{I <: Interpolation}(interp::I, val::Float64)
  if val < interp.x_vals[1]
    return 1
  elseif val >= interp.x_vals[end]
    # return interp.x_vals[end] - interp.x_vals[1] - 2
    return length(interp.x_vals)
  else
    return findfirst(interp.x_vals .> val) - 1 # need to look at this
  end
end

value(interp::LogInterpolation, val::Float64) = exp(value(interp.interpolator, val))

function value(interp::LinearInterpolation, val::Float64)
  i = locate(interp, val)
  # println("I is ", i)
  # println("Val is ", val)
  # println("Y vals ", interp.y_vals)
  # println("X vals ", interp.x_vals)
  # println("S vals ", interp.s)
  return interp.y_vals[i] + (val - interp.x_vals[i]) * interp.s[i]
end
