struct IntegralEngine{B <: AbstractBlackScholesProcess} <: PricingEngine
  process::B
end

struct Integrand{P <: StrikedTypePayoff} <: Function
  payoff::P
  s0::Float64
  drift::Float64
  variance::Float64
end

function (integrand::Integrand)(x::Float64)
  temp = integrand.s0 * exp(x)
  result = integrand.payoff(temp)
  return result * exp(-(x - integrand.drift) * (x - integrand.drift) / (2.0 * integrand.variance))
end

function _calculate!(pe::IntegralEngine, opt::EuropeanOption)
  ex = opt.exercise
  payoff = opt.payoff
  variance = black_variance(pe.process.blackVolatility, ex.dates[end], payoff.strike)

  dividendDiscount = discount(pe.process.dividendYield, ex.dates[end])
  riskFreeDiscount = discount(pe.process.riskFreeRate, ex.dates[end])

  drift_ = log(dividendDiscount / riskFreeDiscount) - 0.5 * variance

  f = Integrand(payoff, state_variable(pe.process).value, drift_, variance)

  integrator = SegmentIntegral(5000)
  infinity = 10.0 * sqrt(variance)

  opt.results.value = discount(pe.process.riskFreeRate, ex.dates[end]) / sqrt(2.0 * pi * variance) *
                      integrator(f, drift_ - infinity, drift_ + infinity)

  return opt
end
