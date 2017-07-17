struct Put <: OptionType end
struct Call <: OptionType end

value(::Put) = -1
value(::Call) = 1

struct PlainVanillaPayoff{OT <: OptionType} <: StrikedTypePayoff
  optionType::OT
  strike::Float64
end

(payoff::PlainVanillaPayoff)(price::Float64) = _get_payoff(payoff.optionType, payoff, price)

_get_payoff(::Call, payoff::PlainVanillaPayoff, price::Float64) = max(price - payoff.strike, 0.0)
_get_payoff(::Put, payoff::PlainVanillaPayoff, price::Float64) = max(payoff.strike - price, 0.0)

struct ForwardTypePayoff{P <: PositionType} <: AbstractPayoff
  position::P
  strike::Float64
end
