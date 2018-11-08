struct BlackCallableFixedRateBondEngine{V <: CallableBondVolatilityStructure} <: PricingEngine
  volatility::V
end

function BlackCallableFixedRateBondEngine(fwdYieldVol::Quote)
  vol = CallableBondConstantVol(0, NullCalendar(), fwdYieldVol, Actual365())

  return BlackCallableFixedRateBondEngine{typeof(vol)}(vol)
end
