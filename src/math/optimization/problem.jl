using QuantLib

mutable struct Problem{F <: CostFunction, C <: Constraint, T}
  costFunction::F
  constraint::C
  initialValue::Vector{T}
  currentValue::Vector{T}
  functionValue::Float64
  squaredNorm::Float64
  functionEvaluation::Int
  gradientEvaluation::Int
end

function Problem(costFunction::F, constraint::C, initialValue::Vector{T}) where {F <: CostFunction, C <: Constraint, T}
  currentValue = initialValue
  functionValue = 0.0
  squaredNorm = 0.0
  functionEvaluation = 0
  gradientEvaluation = 0

  return Problem{F, C, T}(costFunction, constraint, initialValue, currentValue, functionValue, squaredNorm, functionEvaluation, gradientEvaluation)
end

## Problem methods ##
function reset!(p::Problem)
  p.functionEvaluation = p.gradientEvaluation = 0
  p.functionValue = p.squaredNorm = 0.0

  return p
end

function value!(p::Problem, x::Vector{T}) where {T}
  p.functionEvaluation += 1
  return QuantLib.value(p.costFunction, x)
end

function values!(p::Problem, x::Vector{Float64})
  p.functionEvaluation += 1
  return QuantLib.func_values(p.costFunction, x)
end

function extrapolate!(p::Problem, i_highest::Int, factor::Float64, values::Vector{T}, sum_array::Vector{T},
                    vertices::Vector{Vector{T}}) where {T}
  pTry = zeros(length(sum_array))
  while true
    dimensions = length(values) - 1
    factor1 = (1.0 - factor) / dimensions
    factor2 = factor1 - factor
    @inbounds pTry = sum_array * factor1 - vertices[i_highest] * factor2
    factor *= 0.5

    if test(p.constraint, pTry) || abs(factor) <= EPS_VAL
      break
    end
  end

  if abs(factor) <= EPS_VAL
    return values[i_highest]
  end
  factor *= 2
  vTry = value!(p, pTry) # TODO check this
  if vTry < values[i_highest]
    @inbounds values[i_highest] = vTry
    @inbounds sum_array[:] = sum_array + pTry - vertices[i_highest]
    @inbounds vertices[i_highest] = pTry
  end
  # if vTry == 0.02226748007400235
  #   println("sum_array: $sum_array")
  #   println("values: $values")
  #   println("pTry: $pTry")
  #   println("vertices: $vertices")
  #   error("Ooops")
  # end
  return vTry
end
