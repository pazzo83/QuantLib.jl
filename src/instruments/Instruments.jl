## TOP LEVEL CALCULATION METHODS - KEEPING TRACK OF CALCULATION STATE ##
function update_pricing_engine{I <: Instrument, P <: PricingEngine}(inst::I, pe::P)
  T = get_pricing_engine_type(inst)
  if typeof(pe) <: T
    inst.pricingEngine = pe
    inst.lazyMixin.calculated = false
  else
    # we have to clone
    inst = clone(inst; pe = pe)
  end

  return inst
end
