type BaroneAdesiWhaleyApproximationEngine{B <: AbstractBlackScholesProcess} <: PricingEngine
  process::B
end

function get_critical_price_params(::Call, payoff::StrikedTypePayoff, n::Float64, m::Float64,
                                  variance::Float64, bT::Float64)
  qu = (-(n - 1.0) + sqrt(((n - 1.0) * (n - 1.0)) + 4.0 * m)) / 2.0
  Su = payoff.strike / (1.0 - 1.0 / qu)
  h = -(bT + 2.0 * sqrt(variance)) * payoff.strike / (Su - payoff.strike)
  Si = payoff.strike + (Su - payoff.strike) * (1.0 - exp(h))

  return qu, Su, h, Si
end

function get_critical_price_params(::Put, payoff::StrikedTypePayoff, n::Float64, m::Float64,
                                  variance::Float64, bT::Float64)
  qu = (-(n - 1.0) - sqrt(((n - 1.0) * (n - 1.0)) + 4.0 * m)) / 2.0
  Su = payoff.strike / (1.0 - 1.0 / qu)
  h = (bT - 2.0 * sqrt(variance)) * payoff.strike / (payoff.strike - Su)
  Si = Su + (payoff.strike - Su) * exp(h)

  return qu, Su, h, Si
end

function evolve_Si(::Call, payoff::StrikedTypePayoff, Si::Float64, n::Float64, K::Float64, dividendDiscount::Float64,
                  riskFreeDiscount::Float64, tolerance::Float64, variance::Float64, temp::Float64, d1::Float64)
  w = Normal()
  Q = (-(n - 1.0) + sqrt(((n - 1.0) * (n - 1.0)) + 4.0 * K)) / 2.0
  LHS = Si - payoff.strike
  RHS = temp + (1.0 - dividendDiscount * cdf(w, d1)) * Si / Q
  bi = dividendDiscount * cdf(w, d1) * (1.0 - 1.0 / Q) + (1.0 - dividendDiscount * distribution_derivative(w, d1) /
        sqrt(variance)) / Q
  while abs(LHS - RHS) / payoff.strike > tolerance
    Si = (payoff.strike + RHS - bi * Si) / (1.0 - bi)
    forwardSi = Si * dividendDiscount / riskFreeDiscount
    d1 = (log(forwardSi / payoff.strike) + 0.5 * variance) / sqrt(variance)
    LHS = Si - payoff.strike
    temp2 = black_formula(payoff.optionType, payoff.strike, forwardSi, sqrt(variance)) * riskFreeDiscount
    RHS = temp2 + (1.0 - dividendDiscount * cdf(w, d1)) * Si / Q
    bi = dividendDiscount * cdf(w, d1) * (1.0 - 1.0 / Q) + (1.0 - dividendDiscount * distribution_derivative(w, d1) /
          sqrt(variance)) / Q
  end

  return Si
end

function evolve_Si(::Put, payoff::StrikedTypePayoff, Si::Float64, n::Float64, K::Float64, dividendDiscount::Float64,
                  riskFreeDiscount::Float64, tolerance::Float64, variance::Float64, temp::Float64, d1::Float64)
  w = Normal()
  Q = (-(n - 1.0) - sqrt(((n - 1.0) * (n - 1.0)) + 4.0 * K)) / 2.0
  LHS = payoff.strike - Si
  RHS = temp - (1.0 - dividendDiscount * cdf(w, -d1)) * Si / Q
  bi = -dividendDiscount * cdf(w, -d1) * (1.0 - 1.0 / Q) - (1.0 + dividendDiscount * distribution_derivative(w, -d1) /
        sqrt(variance)) / Q
  while abs(LHS - RHS) / payoff.strike > tolerance
    Si = (payoff.strike - RHS + bi * Si) / (1.0 + bi)
    forwardSi = Si * dividendDiscount / riskFreeDiscount
    d1 = (log(forwardSi / payoff.strike) + 0.5 * variance) / sqrt(variance)
    LHS = payoff.strike - Si
    temp2 = black_formula(payoff.optionType, payoff.strike, forwardSi, sqrt(variance)) * riskFreeDiscount
    RHS = temp2 - (1.0 - dividendDiscount * cdf(w, -d1)) * Si / Q
    bi = -dividendDiscount * cdf(w, -d1) * (1.0 - 1.0 / Q) - (1.0 + dividendDiscount * distribution_derivative(w, -d1) /
          sqrt(variance)) / Q
  end

  return Si
end

function critical_price(payoff::StrikedTypePayoff, riskFreeDiscount::Float64, dividendDiscount::Float64,
                        variance::Float64, tolerance::Float64)
  # calculation of seed values, Si
  bT = log(dividendDiscount / riskFreeDiscount)
  n = 2.0 * bT / variance
  m = -2.0 * log(riskFreeDiscount) / variance

  qu, Su, h, Si = get_critical_price_params(payoff.optionType, payoff, n, m, variance, bT)

  # Newton - Raphson algo for finding critical price Si
  forwardSi = Si * dividendDiscount / riskFreeDiscount
  d1 = (log(forwardSi / payoff.strike) + 0.5 * variance) / sqrt(variance)
  K = ~is_close(riskFreeDiscount, 1.0, 1000) ? -2.0 * log(riskFreeDiscount) / (variance * (1.0 - riskFreeDiscount)) :
                                              2.0 / variance

  temp = black_formula(payoff.optionType, payoff.strike, forwardSi, sqrt(variance)) * riskFreeDiscount

  Si = evolve_Si(payoff.optionType, payoff, Si, n, K, dividendDiscount, riskFreeDiscount, tolerance,
                variance, temp, d1)

  return Si
end


function early_exericse_calc(::Call, payoff::StrikedTypePayoff, n::Float64, K::Float64, Sk::Float64,
                            dividendDiscount::Float64, d1::Float64, spot::Float64, black_val::Float64)
  w = Normal()
  Q = (-(n - 1.0) + sqrt(((n - 1.0) * (n - 1.0)) + 4.0 * K)) / 2.0
  a = (Sk / Q) * (1.0 - dividendDiscount * cdf(w, d1))
  if spot < Sk
    return black_val + a * ^((spot / Sk), Q)
  else
    return spot - payoff.strike
  end
end

function early_exericse_calc(::Put, payoff::StrikedTypePayoff, n::Float64, K::Float64, Sk::Float64,
                            dividendDiscount::Float64, d1::Float64, spot::Float64, black_val::Float64)
  w = Normal()
  Q = (-(n - 1.0) - sqrt(((n - 1.0) * (n - 1.0)) + 4.0 * K)) / 2.0
  a = -(Sk / Q) * (1.0 - dividendDiscount * cdf(w, -d1))
  if spot > Sk
    return black_val + a * ^((spot / Sk), Q)
  else
    return payoff.strike - spot
  end
end

function _calculate!(pe::BaroneAdesiWhaleyApproximationEngine, opt::AmericanOption)
  # stuff
  ex = opt.exercise
  payoff = opt.payoff

  variance = black_variance(pe.process.blackVolatility, ex.dates[end], payoff.strike)

  dividendDiscount = discount(pe.process.dividendYield, ex.dates[end])
  riskFreeDiscount = discount(pe.process.riskFreeRate, ex.dates[end])

  spot = state_variable(pe.process).value
  spot > 0.0 || error("negative or null underlying given")
  forwardPrice = spot *  dividendDiscount / riskFreeDiscount

  black = BlackCalculator(payoff, forwardPrice, sqrt(variance), riskFreeDiscount)

  if dividendDiscount > 1.0 && isa(payoff.optionType, Call)
    # early exercise never optimal
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
  else
    # early exercise can be optimal
    tolerance = 1e-6
    Sk = critical_price(payoff, riskFreeDiscount, dividendDiscount, variance, tolerance)
    forwardSk = Sk * dividendDiscount / riskFreeDiscount
    d1 = (log(forwardSk / payoff.strike) + 0.5 * variance) / sqrt(variance)
    n = 2.0 * log(dividendDiscount / riskFreeDiscount)/variance
    K = ~is_close(riskFreeDiscount, 1.0, 1000) ? -2.0 * log(riskFreeDiscount) / (variance * (1.0 - riskFreeDiscount)) :
                                                2.0 / variance
    opt.results.value = early_exericse_calc(payoff.optionType, payoff, n, K, Sk, dividendDiscount, d1,
                          spot, value(black))
  end
end
