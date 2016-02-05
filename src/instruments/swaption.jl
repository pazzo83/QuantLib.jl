type SettlementPhysical <: SettlementType end
type SettlementCash <: SettlementType end

type SwaptionResults{S <: AbstractString}
  value::Float64
  additionalResults::Dict{S, Float64}
end

SwaptionResults() = SwaptionResults(0.0, Dict("spreadCorrection" => 0.0, "strike" => 0.0, "atmForward" => 0.0, "annuity" => 0.0, "swapLength" => 0.0, "stdDev" => 0.0, "vega" => 0.0))

function reset!(results::SwaptionResults)
  results.value = 0.0
  for k in keys(results.additionalResults)
    results.additionalResults[k] = 0.0
  end

  return results
end

type Swaption{E <: Exercise, S <: SettlementType, P <: PricingEngine} <: Option
  lazyMixin::LazyMixin
  swap::VanillaSwap
  exercise::E
  delivery::S
  results::SwaptionResults
  pricingEngine::P

  # function call{E, S}(::Type{Swaption}, lz::LazyMixin, swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults)
  #   new{E, S, PricingEngine}(lz, swap, exercise, delivery, results)
  # end

  function call{E, S, P}(::Type{Swaption}, lz::LazyMixin, swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults, pricingEngine::P = NullSwaptionEngine())
    new{E, S, P}(lz, swap, exercise, delivery, results, pricingEngine)
  end
end

Swaption{E <: Exercise}(swap::VanillaSwap, exercise::E) = Swaption(LazyMixin(), swap, exercise, SettlementPhysical(), SwaptionResults())

function perform_calculations!(swaption::Swaption)
  reset!(swaption.results)
  _calculate!(swaption.pricingEngine, swaption)

  return swaption
end

function npv(swaption::Swaption)
  calculate!(swaption)

  return swaption.results.value
end

function clone(swaption::Swaption; pe::PricingEngine = swaption.pricingEngine)
  lazyMixin, res = pe == swaption.pricingEngine ? (swaption.lazyMixin, swaption.results) : (LazyMixin(), SwaptionResults())

  newSwaption = Swaption(lazyMixin, swaption.swap, swaption.exercise, swaption.delivery, res, pe)

  return newSwaption
end

get_pricing_engine_type{E, S, P}(::Swaption{E, S, P}) = P
