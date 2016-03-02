function black_formula(optType::OptionType, strike::Float64, forward::Float64, stdDev::Float64,
                      discount::Float64 = 1.0, displacement::Float64 = 0.0)
  # TODO check parameters
  if stdDev == 0.0
    return max((forward - strike) * value(optType), 0.0) * discount
  end

  forward = forward + displacement
  strike = strike + displacement

  if strike == 0.0
    return isa(optType, Call) ? forward * discount : 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev
  d2 = d1 - stdDev
  phi = Normal()
  nd1 = phi(value(optionType) * d1)
  nd2 = phi(value(optionType) * d2)
  result = discount * optionType * (forward * nd1 - strike * nd2)
  result >= 0.0 || error("negative value")

  return result
end

function black_scholes_theta(process::AbstractBlackScholesProcess, val::Float64, delta_::Float64, gamma_::Float64)
  u = state_variable(process).value
  r = zero_rate(process.riskFreeRate, 0.0, ContinuousCompounding()).rate
  q = zero_rate(process.dividendYield, 0.0, ContinuousCompounding()).rate
  v = local_vol(process.localVolatility, 0.0, u)

  return r * val - (r - q) * u * delta_ - 0.5 * v * v * u * u * gamma_
end
