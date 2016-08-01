type DiscountingSwapEngine{Y <: YieldTermStructure} <: PricingEngine
  yts::Y
  includeSettlementDateFlows::Bool

  # function call(::Type{DiscountingSwapEngine})
  #   new{YieldTermStructure}()
  # end
  #
  # function call{Y}(::Type{DiscountingSwapEngine}, yts::Y, includeSettlementDateFlows::Bool = true)
  #   new{Y}(yts, includeSettlementDateFlows)
  # end
  DiscountingSwapEngine(yts, includeSettlementDateFlows) = new(yts, includeSettlementDateFlows)
end

DiscountingSwapEngine{Y <: YieldTermStructure}(yts::Y, includeSettlementDateFlows::Bool = true) = DiscountingSwapEngine{Y}(yts, includeSettlementDateFlows)

DiscountingSwapEngine() = DiscountingSwapEngine{NullYieldTermStructure}(NullYieldTermStructure(), true)

function _calculate!{S <: Swap}(pe::DiscountingSwapEngine, swap::S)
  # stuff
  # println("NEW ONE=============================================================================")
  # if swap.rate.value > 0.0323
  #   error("DUIE")
  # end
  swap.results.value = 0.0
  yts = pe.yts

  ref_date = reference_date(yts)
  swap.results.npvDateDiscount = discount(yts, ref_date)
  # for (i, leg) in enumerate(swap.legs)
  for i = 1:length(swap.legs)
    leg = swap.legs[i]
    swap.results.legNPV[i], swap.results.legBPS[i] = npvbps(leg, yts, ref_date, ref_date, pe.includeSettlementDateFlows)
    swap.results.legNPV[i] *= swap.payer[i]
    swap.results.legBPS[i] *= swap.payer[i]

    d1 = accrual_start_date(leg.coupons[1])
    if d1 >= ref_date
      swap.results.startDiscounts[i] = discount(yts, d1)
    end

    d2 = accrual_end_date(leg.coupons[end])
    if (d2 >= ref_date)
      swap.results.endDiscounts[i] = discount(yts, d2)
    end

    swap.results.value += swap.results.legNPV[i]
  end

  if swap.results.legBPS[1] != 0.0
    # println("fixedRate: $(swap.fixedRate)")
    # println("NPV: $(swap.results.value)")
    # println("legBPS: $(swap.results.legBPS[1])")
    swap.results.fairRate = swap.fixedRate - swap.results.value / (swap.results.legBPS[1] / basisPoint)
  end

  return swap
end

clone(pe::DiscountingSwapEngine, ts::TermStructure) = DiscountingSwapEngine(ts)
