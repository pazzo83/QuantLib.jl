struct AnalyticEuropeanEngine{P <: AbstractBlackScholesProcess} <: PricingEngine
  process::P
end

function _calculate!(pe::AnalyticEuropeanEngine, opt::EuropeanOption)
  #TODO check non striked payoff
  payoff = opt.payoff
  variance = black_variance(pe.process.blackVolatility, opt.exercise.dates[end], payoff.strike)

  dividendDiscount = discount(pe.process.dividendYield, opt.exercise.dates[end])

  riskFreeDiscount = discount(pe.process.riskFreeRate, opt.exercise.dates[end])

  spot = state_variable(pe.process).value
  forwardPrice = spot *  dividendDiscount / riskFreeDiscount

  black = BlackCalculator(payoff, forwardPrice, sqrt(variance), riskFreeDiscount)

  opt.results.value = value(black)
  opt.results.delta = delta(black, spot)
  opt.results.deltaForward = delta_forward(black)
  opt.results.elasticity = elasticity(black, spot)
  opt.results.gamma = gamma(black, spot)

  rfdc = pe.process.riskFreeRate.dc
  divdc = pe.process.dividendYield.dc
  voldc = pe.process.blackVolatility.dc

  t = year_fraction(rfdc, reference_date(pe.process.riskFreeRate), opt.exercise.dates[end])
  opt.results.rho = rho(black, t)

  t = year_fraction(divdc, reference_date(pe.process.dividendYield), opt.exercise.dates[end])
  opt.results.dividendRho = dividend_rho(black, t)

  t = year_fraction(voldc, reference_date(pe.process.blackVolatility), opt.exercise.dates[end])
  opt.results.vega = vega(black, t)

  opt.results.theta = theta(black, spot, t)
  opt.results.thetaPerDay = theta_per_day(black, spot, t)

  opt.results.strikeSensitivity = strike_sensitivity(black)
  opt.results.itmCashProbability = itm_cash_probability(black)

  return opt
end
