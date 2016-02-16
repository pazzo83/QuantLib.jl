using JQuantLib

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

abstract CostFunction

abstract Constraint

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

abstract OptimizationMethod

type Simplex <: OptimizationMethod
  lambda::Float64
end

type LevenbergMarquardt <: OptimizationMethod
  epsfcn::Float64
  xtol::Float64
  gtol::Float64
  useCostFunctionsJacobin::Bool
end

LevenbergMarquardt() = LevenbergMarquardt(EPS_VAL, 1.0e-8, 1.0e-8, false)

type Problem{F <: CostFunction, C <: Constraint, T, I <: Integer}
  costFunction::F
  constraint::C
  initialValue::Vector{T}
  currentValue::Vector{T}
  functionValue::Float64
  squaredNorm::Float64
  functionEvaluation::I
  gradientEvaluation::I
end

function Problem{F <: CostFunction, C <: Constraint, T}(costFunction::F, constraint::C, initialValue::Vector{T})
  currentValue = initialValue
  functionValue = 0.0
  squaredNorm = 0.0
  functionEvaluation = 0
  gradientEvaluation = 0

  return Problem(costFunction, constraint, initialValue, currentValue, functionValue, squaredNorm, functionEvaluation, gradientEvaluation)
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



## Levenberg Marquardt Methods ##
function minimize!(lm::LevenbergMarquardt, p::Problem, endCriteria::EndCriteria)
  reset!(p)
  x = p.currentValue
  initCostValues = values!(p, x)

  m = length(initCostValues)
  n = length(x)
  if lm.useCostFunctionsJacobin
    initJacobin = zeros(m, n)
    get_jacobin!(p.costFunction, initJacobin, x)
  end
  fvec = zeros(m)
  diag_ = zeros(n)
  mode = 1
  factor_ = 1.0
  nprint = 0
  info_ = 0
  nfev = 0
  fjac = zeros(m, n)
  ldfjac = m
  ipvt = ones(Int, n)
  qtf = zeros(n)
  wa1 = zeros(n)
  wa2 = zeros(n)
  wa3 = zeros(n)
  wa4 = zeros(m)

  function fcn!{I <: Integer}(::I, n::I, _x::Vector{Float64}, _fvec::Vector{Float64})
    xt = _x[1:n]

    if test(p.constraint, xt)
      tmp = values!(p, xt)
      _fvec[1:length(tmp)] = tmp
    else
      _fvec = copy(initCostValues)
    end

    return _fvec
  end

  # function jacFcn!{I <: Integer}(m::I, n::I, x::Float64, fjac::Float64, ::I)
  #   xt = fill(x + n, n)
  # end

  # function fcn2(x::Vector{Float64})
  #   if test(p.constraint, x)
  #     fvec = values!(p, x)
  #   else
  #     fvec = copy(initCostValues)
  #   end
  #
  #   return fvec
  # end

  # jacFcn! = Calculus.jacobian(fcn2)

  # TODO check requirements, see levenbergmarquardt.cpp
  lmdif!(m, n, x, fvec, endCriteria.functionEpsilon, lm.xtol, lm.gtol, endCriteria.maxIterations, lm.epsfcn, diag_, mode, factor_, nprint, info_, nfev, fjac, ldfjac,
      ipvt, qtf, wa1, wa2, wa3, wa4, fcn!)

  # lmdif2!(n, m, x, mode, factor_, info_, lm.epsfcn, endCriteria.functionEpsilon, lm.xtol, lm.gtol, endCriteria.maxIterations, fcn!)
  # println(x)
  # res = lmfit(fcn2, xx, endCriteria.maxIterations)

  # current value is already set
  p.functionValue = JQuantLib.value(p.costFunction, x)

  return p
end

## SIMPLEX METHODS ##
function compute_simplex_size{T}(vert::Vector{Vector{T}})
  center = sum(vert)
  multiply_array_by_self!(center, 1 / length(vert))

  result = 0.0
  @inbounds for i in vert
    tmp = i - center
    result += sqrt(dot(tmp, tmp))
  end
  return result / length(vert)
end

function minimize!(simplex::Simplex, p::Problem, end_criteria::EndCriteria)
  xtol = end_criteria.rootEpsilon
  max_stationary_state_iterations = end_criteria.maxStationaryStateIterations
  reset!(p)
  x = p.currentValue
  iter_num = 0

  # initialize the vertices of the simplex
  end_condition = false
  n = length(x)
  vertices = Vector{Vector{Float64}}(n + 1)
  fill!(vertices, x)
  direction = zeros(n)
  for i = 1:n
    @inbounds direction[i] = 1.0
    @inbounds vertices[i + 1] = update(p.constraint, vertices[i + 1], direction, simplex.lambda)
    # reset direction
    @inbounds direction[i] = 0.0
  end
  # initialize function values at the vertices of the simplex
  values = zeros(n + 1)
  # values = Vector{DD}(n + 1)
  for i=1:n+1
    @inbounds values[i] = value!(p, vertices[i])
  end

  sum_array = zeros(n + 1)
  # sum_array = Vector{DD}(n + 1)
  # loop through looking for the minimum
  while !end_condition
    # sum_array = zeros(n)
    sum_array = sum(vertices)
    # determine the best (i_lowest) and worst (i_highest) and
    # 2nd worst (i_next_highest) vertices
    i_lowest = 1
    if values[1] < values[2]
      i_highest = 2
      i_next_highest = 1
    else
      i_highest = 1
      i_next_highest = 2
    end

    # we might be able to just do a sort here
    for i=2:n+1
      @inbounds if values[i] > values[i_highest]
        i_next_highest = i_highest
        i_highest = i
      else
        @inbounds if values[i] > values[i_next_highest] && i != i_highest
          i_next_highest = i
        end
      end

      @inbounds if values[i] < values[i_lowest]
        i_lowest = i
      end
    end

    simplex_size = compute_simplex_size(vertices)
    iter_num += 1
    if simplex_size < xtol || iter_num >= end_criteria.maxIterations
      ## TODO implement reason for exiting (EndResult:Type)
      @inbounds x = vertices[i_lowest]
      @inbounds low = values[i_lowest]
      p.functionValue = low
      p.currentValue = x

      return p
    end
    # if end criteria not met, continue
    factor = -1.0
    vTry = extrapolate!(p, i_highest, factor, values, sum_array, vertices)
    @inbounds if vTry <= values[i_lowest] && factor == -1.0
      factor = 2.0
      extrapolate!(p, i_highest, factor, values, sum_array, vertices)
    elseif abs(factor) > EPS_VAL
      @inbounds if vTry >= values[i_next_highest]
        @inbounds vSave = values[i_highest]
        factor = 0.5
        vTry = extrapolate!(p, i_highest, factor, values, sum_array, vertices)
        if vTry >= vSave && abs(factor) > EPS_VAL
          # println("vTry: $vTry")
          # println("vSave: $vSave")
          # println("i_highest: $i_highest")
          # println("i_next_highest: $i_next_highest")
          # println("ending: $iter_num")
          # error("FULL STOP")
          for i = 1:n + 1
            if i != i_lowest
              @inbounds vertices[i] = 0.5 * (vertices[i] + vertices[i_lowest])
              @inbounds values[i] = value!(p, vertices[i])
            end
          end
        end
      end
    end

    if abs(factor) <= EPS_VAL
      x = vertices[i_lowest]
      low = values[i_lowest]
      p.functionValue = low
      p.currentValue = x

      return p
    end
  end
end

centroid(p::Matrix) = reshape(mean(p, 2), size(p, 1))

function dominates(x::Vector, y::Vector)
  for i in 1:length(x)
    @inbounds if x[i] <= y[i]
      return false
    end
  end
  return true
end

function dominates(x::Real, y::Vector)
  for i in 1:length(y)
    @inbounds if x <= y[i]
      return false
    end
  end
  return true
end

nmobjective(y::Vector, m::Integer, n::Integer) = sqrt(var(y) * (m / n))

function minimize_2!(simplex::Simplex, prob::Problem, end_criteria::EndCriteria)
  xtol = ftol = end_criteria.rootEpsilon
  max_stationary_state_iterations = end_criteria.maxStationaryStateIterations
  reset!(prob)
  initial_x = prob.currentValue
  iter_num = 0
  m = length(initial_x)
  initial_step = ones(length(initial_x))
  a = 1.0
  g = 2.0
  b = 0.5

  n = m + 1
  p = repmat(initial_x, 1, n)
  for i = 1:m
    @inbounds p[i, i] += initial_step[i]
  end

  # Maintain a record of the value of f() at n points
  y = Array(Float64, n)
  for i = 1:n
    @inbounds y[i] = value!(prob, p[:, i])
  end

  f_x_previous, f_x = NaN, nmobjective(y, m, n)

  # Count iterations
  iteration = 0

  # Cache p_bar, y_bar, p_star, and p_star_star
  p_bar = Array(Float64, m)
  y_bar = Array(Float64, m)
  p_star = Array(Float64, m)
  p_star_star = Array(Float64, m)

  # Iterate until convergence or exhaustion
  f_converged = false

  while !f_converged && iteration < end_criteria.maxIterations
    # augment iter coutner
    iteration += 1

    # Find p_l and p_h, the min and max values of f() among p
    y_l, l = findmin(y)
    @inbounds p_l = p[:, l]
    y_h, h = findmax(y)
    @inbounds p_h = p[:, h]

    # Compute the centroid of the non-maximal points
    # also cache function values of all non-maximal points
    fill!(p_bar, 0.0)
    tmpindex = 0
    for i = 1:n
      if i != h
        tmpindex += 1
        @inbounds y_bar[tmpindex] = y[i]
        @inbounds p_bar[:] += p[:, i]
      end
    end

    # for j = 1:m
    #   @inbounds p_bar[j] /= m
    # end
    divide_array_by_self!(p_bar, m)

    # Compute a reflection
    for j = 1:m
      @inbounds p_star[j] = (1 + a) * p_bar[j] - a * p_h[j]
    end
    y_star = value!(prob, p_star)

    if y_star < y_l
      # Compute an expansion
      for j = 1:m
        @inbounds p_star_star[j] = g * p_star[j] + (1 - g) * p_bar[j]
      end
      y_star_star = value!(prob, p_star_star)

      if y_star_star < y_l
        @inbounds p_h[:] = p_star_star
        @inbounds p[:, h] = p_star_star
        @inbounds y[h] = y_star_star
      else
        p_h = p_star
        @inbounds p[:, h] = p_star
        @inbounds y[h] = y_star
      end
    else
      if dominates(y_star, y_bar)
        if y_star < y_h
          @inbounds p_h[:] = p_star
          @inbounds p[:, h] = p_h
          @inbounds y[h] = y_star
        end

        # Compute a contraction
        for j = 1:m
          @inbounds p_star_star[j] = b * p_h[j] + (1 - b) * p_bar[j]
        end
        y_star_star = value!(prob, p_star_star)

        if y_star_star > y_h
          for i = 1:n
            for j = 1:m
              @inbounds p[j, i] = (p[j, i] + p_l[j]) / 2.0
            end
            @inbounds y[i] = value!(prob, p[:, i])
          end
        else
          @inbounds p_h[:] = p_star_star
          @inbounds p[:, h] = p_h
          @inbounds y[h] = y_star_star
        end
      else
        @inbounds p_h[:] = p_star
        @inbounds p[:, h] = p_h
        @inbounds y[h] = y_star
      end
    end

    f_x_previous, f_x = f_x, nmobjective(y, m, n)

    if f_x <= ftol
      f_converged = true
    end
  end

  minimum = centroid(p)
  prob.functionValue = Float64(value!(prob, minimum))
  prob.currentValue = minimum

  return prob
end

## Problem methods ##
function reset!(p::Problem)
  p.functionEvaluation = p.gradientEvaluation = 0
  p.functionValue = p.squaredNorm = 0.0

  return p
end

function value!{T}(p::Problem, x::Vector{T})
  p.functionEvaluation += 1
  return JQuantLib.value(p.costFunction, x)
end

function values!(p::Problem, x::Vector{Float64})
  p.functionEvaluation += 1
  return JQuantLib.func_values(p.costFunction, x)
end

function extrapolate!{I <: Integer, T}(p::Problem, i_highest::I, factor::Float64, values::Vector{T}, sum_array::Vector{T},
                    vertices::Vector{Vector{T}})
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
