type G2SwaptionEngine{AffineModelType, Y <: YieldTermStructure, I <: Integer} <: PricingEngine
  model::G2{AffineModelType, Y}
  range::Float64
  intervals::I

  call{Y, I}(::Type{G2SwaptionEngine}, model::G2{AffineModelType, Y}, range::Float64, intervals::I) = new{AffineModelType, Y, I}(model, range, intervals)
end

# methods #
function _calculate!(pe::G2SwaptionEngine, swaption::Swaption)
  swap = swaption.swap

  # overriding pricing engine
  _calculate!(DiscountingSwapEngine(pe.model.ts), swap)
  swap.lazyMixin.calculated = true

  correction = swap.spread * abs(floating_leg_BPS(swap) / fixed_leg_BPS(swap))
  fixedRate = swap.fixedRate - correction
  swaption.results.value = gen_swaption(pe.model, swaption, fixedRate, pe.range, pe.intervals)

  return swaption
end
