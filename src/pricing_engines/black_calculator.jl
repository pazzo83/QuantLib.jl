mutable struct BlackCalculator{S <: StrikedTypePayoff}
  payoff::S
  strike::Float64
  forward::Float64
  stdDev::Float64
  discount::Float64
  variance::Float64
  d1::Float64
  d2::Float64
  alpha::Float64
  beta::Float64
  DalphaDd1::Float64
  DbetaDd2::Float64
  n_d1::Float64
  cum_d1::Float64
  n_d2::Float64
  cum_d2::Float64
  x::Float64
  DxDs::Float64
  DxDstrike::Float64
end

function BlackCalculator(p::S, fwd::Float64, stdDev::Float64, disc::Float64) where {S <: StrikedTypePayoff}
  calc = BlackCalculator{S}(p, p.strike, fwd, stdDev, disc, stdDev * stdDev, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
  initialize!(calc, p)

  return calc
end

function gen_addl_black_vars!(calc::BlackCalculator, optType::Call)
  calc.alpha = calc.cum_d1
  calc.DalphaDd1 = calc.n_d1
  calc.beta = -calc.cum_d2
  calc.DbetaDd2 = -calc.n_d2

  return calc
end

function gen_addl_black_vars!(calc::BlackCalculator, optType::Put)
  calc.alpha = -1.0 + calc.cum_d1
  calc.DalphaDd1 = calc.n_d1
  calc.beta = 1.0 - calc.cum_d2
  calc.DbetaDd2 = -calc.n_d2

  return calc
end

function initialize!(calc::BlackCalculator, p::StrikedTypePayoff)
  calc.strike >= 0.0 || error("strike $strike must be non-negative")
  calc.forward > 0.0 || error("forward $forward must be positive")
  calc.stdDev >= 0.0 || error("stdDev $stdDev must be non-negative")
  calc.discount > 0.0 || error("discount $discount must be positive")

  if calc.stdDev >= eps()
    if is_close(calc.strike, 0.0)
      calc.d1 = JQ_MAX
      calc.d2 = JQ_MAX
      calc.cum_d1 = 1.0
      calc.cum_d2 = 1.0
      calc.n_d1 = 0.0
      calc.n_d2 = 0.0
    else
      calc.d1 = log(calc.forward / calc.strike) / calc.stdDev + 0.5 * calc.stdDev
      calc.d2 = calc.d1 - calc.stdDev
      w = Normal() # Distributions
      calc.cum_d1 = cdf(w, calc.d1)
      calc.cum_d2 = cdf(w, calc.d2)
      calc.n_d1 = distribution_derivative(w, calc.d1)
      calc.n_d2 = distribution_derivative(w, calc.d2)
    end
  else
    if is_close(calc.forward, calc.strike)
      calc.d1 = 0.0
      calc.d2 = 0.0
      calc.cum_d1 = 0.5
      calc.cum_d2 = 0.5
      calc.n_d1 = M_SQRT_2 * M_1_SQRTPI
      calc.n_d2 = M_SQRT_2 * M_1_SQRTPI
    elseif calc.forward > calc.strike
      calc.d1 = JQ_MAX
      calc.d2 = JQ_MAX
      calc.cum_d1 = 1.0
      calc.cum_d2 = 1.0
      calc.n_d1 = 0.0
      calc.n_d2 = 0.0
    else
      calc.d1 = JQ_MIN
      calc.d2 = JQ_MIN
      calc.cum_d1 = 0.0
      calc.cum_d2 = 0.0
      calc.n_d1 = 0.0
      calc.n_d2 = 0.0
    end
  end

  calc.x = calc.strike
  calc.DxDstrike = 1.0

  calc.DxDs = 0.0

  gen_addl_black_vars!(calc, p.optionType)

  return calc
end

value(calc::BlackCalculator) = calc.discount * (calc.forward * calc.alpha + calc.x * calc.beta)

itm_cash_probability(calc::BlackCalculator) = calc.cum_d2
itm_asset_probability(calc::BlackCalculator) = calc.cum_d1

function delta(calc::BlackCalculator, spot::Float64)
  spot > 0.0 || error("positive spot value required.  $spot not allowed")

  DfowardDs = calc.forward / spot

  temp = calc.stdDev * spot
  DalphaDs = calc.DalphaDd1 / temp
  DbetaDs = calc.DbetaDd2 / temp
  temp2 = DalphaDs * calc.forward + calc.alpha * DfowardDs + DbetaDs * calc.x + calc.beta * calc.DxDs

  return calc.discount * temp2
end

function delta_forward(calc::BlackCalculator)
  temp = calc.stdDev * calc.forward
  DalphaDforward = calc.DalphaDd1 / temp
  DbetaDforward = calc.DbetaDd2 / temp
  temp2 = DalphaDforward * calc.forward + calc.alpha + DbetaDforward * calc.x

  return calc.discount * temp2
end

function elasticity(calc::BlackCalculator, spot::Float64)
  val = value(calc)
  del = delta(calc, spot)

  if val > eps()
    return del / val * spot
  elseif abs(del) < eps()
    return 0.0
  elseif del > 0.0
    return JQ_MAX
  else
    return JQ_MIN
  end
end

function gamma(calc::BlackCalculator, spot::Float64)
  spot > 0.0 || error("positive spot value required.  $spot not allowed")

  DforwardDs = calc.forward / spot

  temp = calc.stdDev * spot
  DalphaDs = calc.DalphaDd1 / temp
  DbetaDs = calc.DbetaDd2 / temp

  D2alphaDs2 = -DalphaDs / spot * (1 + calc.d1 / calc.stdDev)
  D2betaDs2 = -DbetaDs / spot * (1 + calc.d2 / calc.stdDev)

  temp2 = D2alphaDs2 * calc.forward + 2.0 * DalphaDs * DforwardDs + D2betaDs2 * calc.x + 2.0 * DbetaDs * calc.DxDs

  return calc.discount * temp2
end

function rho(calc::BlackCalculator, mat::Float64)
  mat >= 0.0 || error("negative maturity not allowed")

  DalphaDr = calc.DalphaDd1 / calc.stdDev
  DbetaDr = calc.DbetaDd2 / calc.stdDev

  temp = DalphaDr * calc.forward + calc.alpha * calc.forward + DbetaDr * calc.x

  return mat * (calc.discount * temp - value(calc))
end

function dividend_rho(calc::BlackCalculator, mat::Float64)
  mat >= 0.0 || error("negative maturity not allowed")

  DalphaDq = -calc.DalphaDd1 / calc.stdDev
  DbetaDq = -calc.DbetaDd2 / calc.stdDev

  temp = DalphaDq * calc.forward - calc.alpha * calc.forward + DbetaDq * calc.x

  return mat * calc.discount * temp
end

function vega(calc::BlackCalculator, mat::Float64)
  temp = log(calc.strike / calc.forward) / calc.variance
  DalphaDsigma = calc.DalphaDd1 * (temp + 0.5)
  DbetaDsigma = calc.DbetaDd2 * (temp - 0.5)

  temp2 = DalphaDsigma * calc.forward + DbetaDsigma * calc.x

  return calc.discount * sqrt(mat) * temp2
end

theta_per_day(calc::BlackCalculator, spot::Float64, mat::Float64) = theta(calc, spot, mat) / 365.0

function theta(calc::BlackCalculator, spot::Float64, mat::Float64)
  if mat < 0.0
    println("WARN: Maturity is negative")
    return -1.0
  end

  is_close(mat, 0.0) && return 0.0

  return -( log(calc.discount) * value(calc) + log(calc.forward / spot) * spot * delta(calc, spot)
          + 0.5 * calc.variance * spot * spot * gamma(calc, spot)) / mat
end

function strike_sensitivity(calc::BlackCalculator)
  temp = calc.stdDev * calc.strike
  DalphaDstrike = -calc.DalphaDd1 / temp
  DbetaDstrike = -calc.DbetaDd2 / temp

  temp2 = DalphaDstrike * calc.forward + DbetaDstrike * calc.x + calc.beta * calc.DxDstrike

  return calc.discount * temp2
end
