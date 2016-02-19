using JQuantLib

abstract OptimizationMethod
abstract CostFunction
abstract Constraint

const FINITE_DIFFERENCES_EPSILON = 1e-8

type Projection{I <: Integer}
  actualParameters::Vector{Float64}
  fixedParameters::Vector{Float64}
  fixParams::BitArray
  numberOfFreeParams::I
end

function Projection(parameterValues::Vector{Float64}, fixParams::BitArray)
  # get num of free params
  numFree = 0
  for i in fixParams
    if !i
      numFree += 1
    end
  end

  return Projection(parameterValues, parameterValues, fixParams, numFree)
end

function project(proj::Projection, params::Vector{Float64})
  projectedParams = Vector{Float64}(proj.numberOfFreeParams)

  i = 1
  for j = 1:length(proj.fixParams)
    if !proj.fixParams[j]
      projectedParams[i] = params[j]
      i += 1
    end
  end

  return projectedParams
end

function include_params(proj::Projection, params::Vector{Float64})
  y = copy(proj.fixedParameters)
  i = 1
  for j = 1:length(y)
    if !proj.fixParams[j]
      y[j] = params[i]
      i += 1
    end
  end
  return y
end

type NoConstraint <: Constraint end
type PositiveConstraint <: Constraint end
type BoundaryConstraint <: Constraint
  low::Float64
  high::Float64
end

type ProjectedConstraint{C <: Constraint} <: Constraint
  constraint::C
  projection::Projection
end

type EndCriteria{I <: Integer}
  maxIterations::I
  maxStationaryStateIterations::I
  rootEpsilon::Float64
  functionEpsilon::Float64
  gradientNormEpsilon::Float64
end

test{T}(::NoConstraint, ::Vector{T}) = true

test(c::ProjectedConstraint, x::Vector{Float64}) = test(c.constraint, include_params(c.projection, x))

function test(::PositiveConstraint, x::Vector{Float64})
  for i = 1:length(x)
    if x[i] <= 0.0
      return false
    end
  end

  return true
end

function test(c::BoundaryConstraint, x::Vector{Float64})
  for i = 1:length(x)
    if x[i] < c.low || x[i] > c.high
      return false
    end
  end

  return true
end

function update{C <: Constraint, T}(constraint::C, params::Vector{T}, direction::Vector{Float64}, beta::Float64)
  diff = beta
  new_params = params + diff * direction
  valid = test(constraint, new_params)
  icount = 0
  while !valid
    if (icount > 200)
      error("Can't update parameter vector")
    end

    diff *= 0.5
    icount += 1
    new_params = params + diff * direction
    valid = test(constraint, new_params)
  end

  params += diff * direction
  return params
end

## Cost Function methods ##
function get_jacobin!{C <: CostFunction}(cf::C, jac::Matrix{Float64}, x::Vector{Float64})
  eps_ = FINITE_DIFFERENCES_EPSILON
  xx = zeros(length(x))
  for i = 1:length(x)
    xx[i] += eps_
    fp = JQuantLib.func_values(cf, xx)
    xx[i] -= 2.0 * eps_
    fm = JQuantLib.func_values(cf, xx)
    for j = 1:length(fp)
      jac[j,i] = 0.5 * (fp[j] - fm[j]) / eps_
    end

    xx[i] = x[i]
  end
  return jac
end
