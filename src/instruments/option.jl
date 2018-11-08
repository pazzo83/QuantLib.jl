# type aliases
const EuropeanOption = Option{EuropeanExercise}
const AmericanOption = Option{AmericanExercise}
const BermudanOption = Option{BermudanExercise}

mutable struct OptionResults # Greeks
  delta::Float64
  gamma::Float64
  theta::Float64
  vega::Float64
  rho::Float64
  dividendRho::Float64
  deltaForward::Float64
  elasticity::Float64
  thetaPerDay::Float64
  strikeSensitivity::Float64
  itmCashProbability::Float64
  value::Float64
end

OptionResults() = OptionResults(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

struct VanillaOptionArgs{P <: StrikedTypePayoff, E <: Exercise}
  payoff::P
  exercise::E
end

function reset!(res::OptionResults)
  res.delta = 0.0
  res.gamma = 0.0
  res.theta = 0.0
  res.vega = 0.0
  res.rho = 0.0
  res.dividendRho = 0.0
  res.deltaForward = 0.0
  res.elasticity = 0.0
  res.thetaPerDay = 0.0
  res.strikeSensitivity = 0.0
  res.itmCashProbability = 0.0
  res.value = 0.0

  return res
end

mutable struct VanillaOption{S <: StrikedTypePayoff, E <: Exercise, P <: PricingEngine} <: OneAssetOption{E}
  lazyMixin::LazyMixin
  payoff::S
  exercise::E
  pricingEngine::P
  results::OptionResults
end

VanillaOption(payoff::S, exercise::E, pe::P) where {S <: StrikedTypePayoff, E <: Exercise, P <: PricingEngine} =
              VanillaOption{S, E, P}(LazyMixin(), payoff, exercise, pe, OptionResults())

get_pricing_engine_type(::VanillaOption{S, E, P}) where {S, E, P} = P

function clone(opt::VanillaOption, pe::P = opt.pricingEngine) where {P <: PricingEngine}
  lazyMixin, res = pe == opt.pricingEngine ? (opt.lazyMixin, opt.results) : (LazyMixin(), OptionResults())

  return VanillaOption{typeof(opt.payoff), typeof(opt.exercise), P}(lazyMixin, opt.payoff, opt.exercise, pe, res)
end
