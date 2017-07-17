struct BjerksundStenslandApproximationEngine{B <: AbstractBlackScholesProcess} <: PricingEngine
  process::B
end

function get_phi(S::Float64, gamma::Float64, H::Float64, I::Float64, rT::Float64, bT::Float64, variance::Float64)
  lambda = (-rT + gamma * bT + 0.5 * gamma * (gamma - 1.0) * variance)
  d = -(log(S / H) + (bT + (gamma - 0.5) * variance)) / sqrt(variance)
  kappa = 2.0 * bT / variance + (2.0 * gamma - 1.0)
  w = Normal()
  return exp(lambda) * (cdf(w, d) - ^((I / S), kappa) * cdf(w, d - 2.0 * log(I / S) / sqrt(variance)))
end

function american_call_approximation(S::Float64, X::Float64, rfD::Float64, dD::Float64, variance::Float64)
  bT = log(dD / rfD)
  rT = log(1.0 / rfD)

  beta = (0.5 - bT / variance) + sqrt(^((bT / variance - 0.5), 2.0) + 2.0 * rT / variance)
  BInfinity = beta / (beta - 1.0) * X
  B0 = max(X, rT / (rT - bT) * X)
  ht = -(bT + 2.0 * sqrt(variance)) * B0 / (BInfinity - B0)

  # investigate what happens to I for db -> 0.0
  I = B0 + (BInfinity - B0) * (-expm1(ht))
  I >= X || error("Bjerksund - Stensland approximation not applicable to this set of params")

  if S >= I
    return S - X
  else
    # investigate what happens to alpha for dD -> 0.0
    return (I - X) * ^(S / I, beta) * (1.0 - get_phi(S, beta, I, I, rT, bT, variance)) +
        S * get_phi(S, 1.0, I, I, rT, bT, variance) -
        S * get_phi(S, 1.0, X, I, rT, bT, variance) -
        X * get_phi(S, 0.0, I, I, rT, bT, variance) +
        X * get_phi(S, 0.0, X, I, rT, bT, variance)
  end
end

function _calculate!(pe::BjerksundStenslandApproximationEngine, opt::AmericanOption)
  ex = opt.exercise
  payoff = opt.payoff

  variance = black_variance(pe.process.blackVolatility, ex.dates[end], payoff.strike)

  dividendDiscount = discount(pe.process.dividendYield, ex.dates[end])
  riskFreeDiscount = discount(pe.process.riskFreeRate, ex.dates[end])

  spot = state_variable(pe.process).value
  spot > 0.0 || error("negative or null underlying given")
  strike = payoff.strike

  if isa(payoff.optionType, Put)
    # use put-call symmetry
    spot, strike = strike, spot
    riskFreeDiscount, dividendDiscount = dividendDiscount, riskFreeDiscount

    payoff = PlainVanillaPayoff(Call(), strike)
  end

  if dividendDiscount >= 1.0
    # early exercise is never optimal - use Black formula
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
  else
    # early exericse can be optimal - use approximation
    opt.results.value = american_call_approximation(spot, strike, riskFreeDiscount, dividendDiscount, variance)
  end

  return opt
end
