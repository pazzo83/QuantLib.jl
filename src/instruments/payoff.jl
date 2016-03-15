type Put <: OptionType end
type Call <: OptionType end

value(::Put) = -1
value(::Call) = 1

type PlainVanillaPayoff{OT <: OptionType} <: StrikedTypePayoff
  optionType::OT
  strike::Float64
end

call(payoff::PlainVanillaPayoff, price::Float64) = _get_payoff(payoff.optionType, payoff, price)

_get_payoff(::Call, payoff::PlainVanillaPayoff, price::Float64) = max(price - payoff.strike, 0.0)
_get_payoff(::Put, payoff::PlainVanillaPayoff, price::Float64) = max(payoff.strike - price, 0.0)

type ForwardTypePayoff{P <: PositionType} <: AbstractPayoff
  position::P
  strike::Float64
end
