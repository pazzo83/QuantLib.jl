type GeneralBlackScholesType <: BlackScholesType end
type BlackScholesMertonType <: BlackScholesType end

type GeneralizedBlackScholesProcess{Y1 <: YieldTermStructure, Y2 <: YieldTermStructure, B <: BlackVolTermStructure, D <: AbstractDiscretization, BST <: BlackScholesType} <: AbstractBlackScholesProcess
  x0::Quote
  riskFreeRate::Y1
  dividendYield::Y2
  blackVolatility::B
  localVolatility::LocalConstantVol
  disc::D
  isStrikeDependent::Bool
  blackScholesType::BST
end

function GeneralizedBlackScholesProcess(x0::Quote, riskFreeRate::YieldTermStructure, dividendYield::YieldTermStructure,
          blackVolatility::BlackConstantVol, disc::AbstractDiscretization = EulerDiscretization())
  localVolatility = LocalConstantVol(blackVolatility.referenceDate, black_vol(blackVolatility, 0.0, x0.value), blackVolatility.dc)

  return GeneralizedBlackScholesProcess(x0, riskFreeRate, dividendYield, blackVolatility, localVolatility, disc, true, GeneralBlackScholesType())
end

# type BlackScholesMertonProcess{Y1 <: YieldTermStructure, Y2 <: YieldTermStructure, B <: BlackVolTermStructure, D <: AbstractDiscretization} <: AbstractBlackScholesProcess
#   x0::Quote
#   riskFreeRate::Y1
#   dividendYield::Y2
#   blackVolatility::B
#   localVolatility::LocalConstantVol
#   disc::D
#   isStrikeDependent::Bool
# end

function BlackScholesMertonProcess(x0::Quote, riskFreeRate::YieldTermStructure, dividendYield::YieldTermStructure,
          blackVolatility::BlackConstantVol, disc::AbstractDiscretization = EulerDiscretization())
  localVolatility = LocalConstantVol(blackVolatility.referenceDate, black_vol(blackVolatility, 0.0, x0.value), blackVolatility.dc)

  return GeneralizedBlackScholesProcess(x0, riskFreeRate, dividendYield, blackVolatility, localVolatility, disc, true, BlackScholesMertonType())
end

function drift(process::GeneralizedBlackScholesProcess, t::Float64, x::Float64)
  sigma = diffusion(process, t, x)
  t1 = t + 0.0001

  return forward_rate(process.riskFreeRate, t, t1, ContinuousCompounding(), NoFrequency()).rate -
          forward_rate(process.dividendYield, t, t1, ContinuousCompounding(), NoFrequency()).rate -
          0.5 * sigma * sigma
end

diffusion(process::GeneralizedBlackScholesProcess, t::Float64, x::Float64) = local_vol(process.localVolatility, t, x)

apply(process::GeneralizedBlackScholesProcess, x0::Float64, dx::Float64) = x0 * exp(dx)

expectation(process::GeneralizedBlackScholesProcess, ::Float64, ::Float64, ::Float64) = error("not implemented")

get_x0(process::AbstractBlackScholesProcess) = process.x0.value

function evolve(process::GeneralizedBlackScholesProcess, t0::Float64, x0::Float64, dt::Float64, dw::Float64)
  if process.isStrikeDependent
    var = black_variance(process.blackVolatility, t0 + dt, 0.01) - black_variance(process.blackVolatility, t0, 0.01)
    drift_ = (forward_rate(process.riskFreeRate, t0, t0 + dt, ContinuousCompounding(), NoFrequency()).rate -
              forward_rate(process.dividendYield, t0, t0 + dt, ContinuousCompounding(), NoFrequency()).rate) *
              dt - 0.5 * var

    return x0 * exp(sqrt(var) * dw + drift_)
  else
    return apply(process, x0, drift(process.disc, process, t0, x0, dt) + std_deviation(process, t0, x0, dt) * dw)
  end
end

state_variable(process::AbstractBlackScholesProcess) = process.x0

get_time(process::AbstractBlackScholesProcess, d::Date) = year_fraction(process.riskFreeRate.dc, reference_date(process.riskFreeRate), d)
