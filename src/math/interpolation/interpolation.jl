abstract Interpolation
abstract Interpolation2D <: Interpolation

abstract DerivativeApprox
abstract BoundaryCondition

type Spline <: DerivativeApprox end
type Lagrange <: BoundaryCondition end

# general interpolation methods #
# update if value passed in
function update!(interp::Interpolation, idx::Int, val::Float64)
  interp.y_vals[idx] = val

  update!(interp, idx)

  return interp
end

# locate x
function locate(interp::Interpolation, val::Float64)
  if val < interp.x_vals[1]
    return 1
  elseif val >= interp.x_vals[end - 1]
    # return interp.x_vals[end] - interp.x_vals[1] - 2
    return length(interp.x_vals) - 1
  else
    # return findfirst(interp.x_vals .> val) - 1 # need to look at this
    return searchsortedlast(interp.x_vals, val)
  end
end
