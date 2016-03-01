type BlackScholesMertonProcess{Y1 <: YieldTermStructure, Y2 <: YieldTermStructure, B <: BlackVolTermStructure, D <: AbstractDiscretization} <: AbstractBlackScholesProcess
  x0::Quote
  riskFreeRate::Y1
  dividendYield::Y2
  blackVolatility::B
  localVolatility::LocalConstantVol
  disc::D
  isStrikeDependent::Bool
end

function BlackScholesMertonProcess(x0::Quote, riskFreeRate::YieldTermStructure, dividendYield::YieldTermStructure,
          blackVolatility::BlackConstantVol, disc::AbstractDiscretization = EulerDiscretization())
  localVolatility = LocalConstantVol(blackVolatility.referenceDate, black_vol(blackVolatility, 0.0, x0.value), blackVolatility.dc)

  return BlackScholesMertonProcess(x0, riskFreeRate, dividendYield, blackVolatility, localVolatility, disc, true)
end

function drift(process::BlackScholesMertonProcess, t::Float64, x::Float64)
  sigma = diffusion(process, t, x)
  t1 = t + 0.0001

  return forward_rate(process.riskFreeRate, t, t1, ContinuousCompounding(), NoFrequency()).rate -
          forward_rate(process.dividendYield, t, t1, ContinuousCompounding(), NoFrequency()).rate -
          0.5 * sigma * sigma
end

diffusion(process::BlackScholesMertonProcess, t::Float64, x::Float64) = local_vol(process.localVolatility, t, x)

apply(process::BlackScholesMertonProcess, x0::Float64, dx::Float64) = x0 * exp(dx)

expectation(process::BlackScholesMertonProcess, ::Float64, ::Float64, ::Float64) = error("not implemented")

function evolve(process::BlackScholesMertonProcess, t0::Float64, x0::Float64, dt::Float64, dw::Float64)
  if process.isStrikeDependent
    var = black_variance(process.blackVolatility, t0 + dt, 0.01) - black_variance(process.blackVolatility, t0, 0.01)
    drift_ = (forward_rate(process.riskFreeRate, t0, t0 + dt, ContinuousCompounding(), NoFrequency()) -
              forward_rate(process.dividendYield, t0, t0 + dt, ContinuousCompounding(), NoFrequency())) *
              dt - 0.5 * var

    return x0 * exp(sqrt(var) * dw + drift)
  else
    return apply(process, x0, drift(process.disc, process, t0, x0, dt) + std_deviation(process, t0, x0, dt) * dw)
  end
end

state_variable(process::AbstractBlackScholesProcess) = process.x0

get_time(process::AbstractBlackScholesProcess, d::Date) = year_fraction(process.riskFreeRate.dc, reference_date(process.riskFreeRate), d)
