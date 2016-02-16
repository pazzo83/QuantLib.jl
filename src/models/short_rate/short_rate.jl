## ShortRateTree methods ##
get_state_prices!(tree::ShortRateTree, i::Int) = get_state_prices!(tree.treeLattice, i)

type PrivateConstraint{P <: Parameter} <: Constraint
  arguments::Vector{P}
end

# type ShortRateModelCommon
#   observers::Vector
# end
#
# ShortRateModelCommon() = ShortRateModelCommon([])

function JQuantLib.Math.test(c::PrivateConstraint, x::Vector{Float64})
  k = 1
  for i = 1:length(c.arguments)
    sz = length(c.arguments[i].data)
    testParams = zeros(sz)
    for j = 1:sz
      testParams[j] = x[k]
      k += 1
    end
    if !test_params(c.arguments[i], testParams)
      return false
    end
  end

  return true
end

type CalibrationFunction{M <: ShortRateModel, C <: CalibrationHelper} <: CostFunction
  model::M
  helpers::Vector{C}
  weights::Vector{Float64}
  projection::Projection
end

type SolvingFunction
  lambda::Vector{Float64}
  Bb::Vector{Float64}
end

function operator(solvFunc::SolvingFunction)
  function _inner(y::Float64)
    value_ = 1.0
    for i = 1:length(solvFunc.lambda)
      value_ -= solvFunc.lambda[i] * exp(-solvFunc.Bb[i] * y)
    end
    return value_
  end

  return _inner
end

function func_values(calibF::CalibrationFunction, params::Vector{Float64})
  set_params!(calibF.model, include_params(calibF.projection, params))
  values = zeros(length(calibF.helpers))
  for i = 1:length(values)
    values[i] = calibration_error(calibF.helpers[i].calibCommon.calibrationErrorType, calibF.helpers[i]) * sqrt(calibF.weights[i])
  end

  return values
end

function value(calibF::CalibrationFunction, params::Vector{Float64})
  set_params!(calibF.model, include_params(calibF.projection, params))
  _value = 0.0
  for i = 1:length(calibF.helpers)
    diff = calibration_error(calibF.helpers[i].calibCommon.calibrationErrorType, calibF.helpers[i])
    _value += diff * diff * calibF.weights[i]
  end

  return sqrt(_value)
end

# accessor methods ##
get_a{M <: ShortRateModel}(m::M) = m.a.data[1]
get_sigma{M <: ShortRateModel}(m::M) = m.sigma.data[1]
get_b{M <: ShortRateModel}(m::M) = m.b.data[1]

check_params_equal(m::ShortRateModel, params::Vector{Float64}) = get_params(m) == params

function set_params!{M <: ShortRateModel}(model::M, params::Vector{Float64})
  if check_params_equal(model, params)
    return model
  end
  paramCount = 1
  args = model.privateConstraint.arguments
  for i = 1:length(args)
    arg = model.privateConstraint.arguments[i]
    for j = 1:length(arg.data)
       arg.data[j] = params[paramCount]
       paramCount += 1
     end
   end
   generate_arguments!(model)

   notify_observers!(model)
   # model.phi = G2FittingParameter(get_a(model), get_sigma(model), get_b(model), get_eta(model), get_rho(model), model.ts)
   return model
 end

## Short Rate Model calibration function #
function calibrate!{M <: ShortRateModel, C <: CalibrationHelper, O <: OptimizationMethod}(model::M, instruments::Vector{C}, method::O, endCriteria::EndCriteria,
                    constraint::Constraint = model.privateConstraint, weights::Vector{Float64} = ones(length(instruments)), fixParams::BitArray = BitArray(0))

  w = length(weights) == 0 ? ones(length(instruments)) : weights
  prms = get_params(model)
  println("model params: ", prms)
  all = falses(length(prms))
  proj = Projection(prms, length(fixParams) > 0 ? fixParams : all)
  calibFunc = CalibrationFunction(model, instruments, w, proj)
  pc = ProjectedConstraint(constraint, proj)
  prob = Problem(calibFunc, pc, project(proj, prms))

  # minimization
  minimize!(method, prob, endCriteria)
  res = prob.currentValue
  set_params!(model, include_params(proj, res))

  return model
end
