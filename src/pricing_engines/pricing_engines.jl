# Pricing Engines module
# module PricingEngines

type DiscountingBondEngine{Y <: YieldTermStructure} <: PricingEngine
  yts::Y

  function call(::Type{DiscountingBondEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingBondEngine}, yts::Y)
    new{Y}(yts)
  end
end

type DiscountingSwapEngine{Y <: YieldTermStructure} <: PricingEngine
  yts::Y

  function call(::Type{DiscountingSwapEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingSwapEngine}, yts::Y)
    new{Y}(yts)
  end
end

function _calculate!{B <: Bond}(pe::DiscountingBondEngine, bond::B)
  yts = pe.yts
  valuation_date = yts.referenceDate
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)
  bond.settlementValue = value
  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return bond
end

function _calculate!{S <: Swap}(pe::DiscountingSwapEngine, swap::S)
  # stuff
  # println("NEW ONE=============================================================================")
  # if swap.rate.value > 0.0323
  #   error("DUIE")
  # end
  swap.results.value = 0.0
  yts = pe.yts

  ref_date = yts.referenceDate

  swap.results.npvDateDiscount = discount(yts, ref_date)

  # for (i, leg) in enumerate(swap.legs)
  for i = 1:length(swap.legs)
    leg = swap.legs[i]
    swap.results.legNPV[i], swap.results.legBPS[i] = npvbps(leg, yts, ref_date, ref_date)
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
    swap.results.fairRate = swap.fixedRate - swap.results.value / (swap.results.legBPS[1] / basisPoint)
  end

  return swap
end
