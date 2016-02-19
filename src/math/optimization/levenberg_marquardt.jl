type LevenbergMarquardt <: OptimizationMethod
  epsfcn::Float64
  xtol::Float64
  gtol::Float64
  useCostFunctionsJacobin::Bool
end

LevenbergMarquardt() = LevenbergMarquardt(EPS_VAL, 1.0e-8, 1.0e-8, false)

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
