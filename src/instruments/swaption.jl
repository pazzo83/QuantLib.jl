struct SettlementPhysical <: SettlementType end
struct SettlementCash <: SettlementType end

mutable struct SwaptionResults
  value::Float64
  additionalResults::Dict{String, Float64}
end

SwaptionResults() = SwaptionResults(0.0, Dict("spreadCorrection" => 0.0, "strike" => 0.0, "atmForward" => 0.0, "annuity" => 0.0, "swapLength" => 0.0, "stdDev" => 0.0,
                                              "vega" => 0.0))

function reset!(results::SwaptionResults)
  results.value = 0.0
  for k in keys(results.additionalResults)
    results.additionalResults[k] = 0.0
  end

  return results
end

mutable struct Swaption{E <: Exercise, S <: SettlementType, P <: PricingEngine, VS <: VanillaSwap} <: Option{E}
  lazyMixin::LazyMixin
  swap::VS
  exercise::E
  delivery::S
  results::SwaptionResults
  pricingEngine::P

  # function call{E, S}(::Type{Swaption}, lz::LazyMixin, swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults)
  #   new{E, S, PricingEngine}(lz, swap, exercise, delivery, results)
  # end

  # function call{E, S, P}(::Type{Swaption}, lz::LazyMixin, swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults, pricingEngine::P = NullSwaptionEngine())
  #   new{E, S, P}(lz, swap, exercise, delivery, results, pricingEngine)
  # end
  # Swaption{E, S}(lz::LazyMixin, swap::VanillaSwap, exercise::E, delivery::S, results::SwaptionResults, pricingEngine::P = NullSwaptionEngine()) =
  #   new{E, S, P}(lz, swap, exercise, delivery, results, pricingEngine)
end

Swaption(swap::VS, exercise::E) where {VS <: VanillaSwap, E <: Exercise} = Swaption{E, SettlementPhysical, NullSwaptionEngine, VS}(LazyMixin(), swap, exercise, SettlementPhysical(), SwaptionResults(), NullSwaptionEngine())

Swaption(swap::VS, exercise::E, pe::P) where {E <: Exercise, P <: PricingEngine, VS <: VanillaSwap} = Swaption{E, SettlementPhysical, P, VS}(LazyMixin(), swap, exercise, SettlementPhysical(), SwaptionResults(), pe)

mutable struct NonstandardSwaption{E <: Exercise, S <: SettlementType, P <: PricingEngine, NS <: NonstandardSwap} <: Option{E}
  lazyMixin::LazyMixin
  swap::NS
  exercise::E
  delivery::S
  results::SwaptionResults
  pricingEngine::P

  # function call{E, S, P}(::Type{NonstandardSwaption}, lz::LazyMixin, swap::NonstandardSwap, exercise::E, delivery::S, results::SwaptionResults, pricingEngine::P = NullSwaptionEngine())
  #   new{E, S, P}(lz, swap, exercise, delivery, results, pricingEngine)
  # end
end

NonstandardSwaption(swap::NS, exercise::E) where {E <: Exercise, NS <: NonstandardSwap} = NonstandardSwaption{E, SettlementPhysical, NullSwaptionEngine, NS}(LazyMixin(), swap, exercise,
                                                                      SettlementPhysical(), SwaptionResults(), NullSwaptionEngine())
NonstandardSwaption(swap::NS, exercise::E, pe::P) where {E <: Exercise, P <: PricingEngine, NS <: NonstandardSwap} =
    NonstandardSwaption{E, SettlementPhysical, P, NS}(LazyMixin(), swap, exercise, SettlementPhysical(), SwaptionResults(), pe)

function perform_calculations!(swaption::Union{Swaption, NonstandardSwaption})
  reset!(swaption.results)
  _calculate!(swaption.pricingEngine, swaption)

  return swaption
end

function npv(swaption::Union{Swaption, NonstandardSwaption})
  calculate!(swaption)

  return swaption.results.value
end

function calibration_basket(swaption::NonstandardSwaption, swaptionEngine::PricingEngine, swapIndex::SwapIndex, swaptionVol::SwaptionVolatilityStructure,
                            basketType::CalibrationBasketType)
  return calibration_basket(swaptionEngine, swaption, swapIndex, swaptionVol, basketType)
end

function clone(swaption::Swaption, pe::PricingEngine = swaption.pricingEngine)
  lazyMixin, res = pe == swaption.pricingEngine ? (swaption.lazyMixin, swaption.results) : (LazyMixin(), SwaptionResults())

  newSwaption = Swaption(lazyMixin, swaption.swap, swaption.exercise, swaption.delivery, res, pe)

  return newSwaption
end

get_pricing_engine_type(::Swaption{E, S, P, NS}) where {E, S, P, NS} = P
