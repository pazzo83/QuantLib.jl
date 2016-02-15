# Pricing Engines module
# module PricingEngines

type NullPricingEngine <: PricingEngine end

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
  includeSettlementDateFlows::Bool

  function call(::Type{DiscountingSwapEngine})
    new{YieldTermStructure}()
  end

  function call{Y}(::Type{DiscountingSwapEngine}, yts::Y, includeSettlementDateFlows::Bool = true)
    new{Y}(yts, includeSettlementDateFlows)
  end
end

type MidPointCdsEngine{P <: AbstractDefaultProbabilityTermStructure, Y <: YieldTermStructure} <: PricingEngine
  probability::P
  recoveryRate::Float64
  discountCurve::Y
end

_calculate!(pe::NullPricingEngine, inst::Instrument) = error("You must set a valid pricing engine")

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

function _calculate!(pe::MidPointCdsEngine, swap::CreditDefaultSwap)
  today_ = settings.evaluation_date
  settlementDate = reference_date(pe.discountCurve)

  # Upfront flow NPV, either we are on the run (no flow) or we are forward start
  upFPV01 = 0.0
  if has_occurred(swap.upfrontPayment, settlementDate)
    # date determining the probability survival so we have to pay the upfront (did not knock out)
    effectiveUpfrontDate = swap.protectionStart > reference_date(pe.probability) ? swap.protectionStart : reference_date(pe.probability)

    upFPV01 = survival_probability(pe.probability, effectiveUpfrontDate) * discount(pe.discountCurve, date(swap.upfrontPayment))
  end

  swap.results.upfrontNPV = upFPV01 * amount(swap.upfrontPayment)
  swap.results.couponLegNPV = 0.0
  swap.results.defaultLegNPV = 0.0

  for i in eachindex(swap.leg.coupons)
    if has_occurred(swap.leg.coupons[i], settlementDate)
      continue
    end

    coup = swap.leg.coupons[i]

    # first calculate NPV of both legs as positive, applying correct sign at the end
    paymentDate = date(coup)
    startDate = accrual_start_date(coup)
    endDate = accrual_end_date(coup)

    if i == 1
      startDate = swap.protectionStart
    end

    effectiveStartDate = (startDate <= today_ && today_ <= endDate) ? today_ : startDate
    defaultDate = effectiveStartDate + Dates.Day(floor(Int(endDate - effectiveStartDate) / 2)) # midpoint

    S = survival_probability(pe.probability, paymentDate)
    P = default_probability(pe.probability, effectiveStartDate, endDate)

    # on one side, we add the fixed rate payments in case of survival
    swap.results.couponLegNPV += S * amount(coup) * discount(pe.discountCurve, paymentDate)

    # possibly including accrual in case of default
    if swap.settlesAccrual
      if swap.paysAtDefaultTime
        swap.results.couponLegNPV += P * accrued_amount(coup, defaultDate) * discount(pe.discountCurve, defaultDate)
      else
        # pays at end
        swap.results.couponLegNPV += P * amount(coup) * discount(pe.discountCurve, paymentDate)
      end
    end

    # on the other side we add the payment in case of default
    claim = amount(swap.claim, defaultDate, swap.notional, pe.recoveryRate)

    if swap.paysAtDefaultTime
      swap.results.defaultLegNPV += P * claim * discount(pe.discountCurve, defaultDate)
    else
      swap.results.defaultLegNPV += P * claim * discount(pe.discountCurve, paymentDate)
    end

  end

  upfrontSign = 1.0
  if isa(swap.side, Seller)
    swap.results.defaultLegNPV *= -1.0
  else
    swap.results.couponLegNPV *= -1.0
    swap.results.upfrontNPV *= -1.0
    upfrontSign = -1.0
  end

  swap.results.value = swap.results.defaultLegNPV + swap.results.couponLegNPV + swap.results.upfrontNPV

  if swap.results.couponLegNPV != 0.0
    swap.results.fairSpread = -swap.results.defaultLegNPV * swap.spread / swap.results.couponLegNPV
  else
    swap.results.fairSpread = -1.0
  end

  upfrontSensitivity = upFPV01 * swap.notional
  if upfrontSensitivity != 0.0
    swap.results.fairUpfront = -upfrontSign * (swap.results.defaultLegNPV + swap.results.couponLegNPV) / upfrontSensitivity
  else
    swap.results.fairUpfront = -1.0
  end

  basisPoint = 1.0e-4
  if swap.spread != 0.0
    swap.results.couponLegBPS = swap.results.couponLegNPV * basisPoint / swap.spread
  else
    swap.results.couponLegBPS = -1.0
  end

  # need to factor in swap upfront
  swap.results.upfrontBPS = -1.0

  return swap
end

# clone methods #
clone(pe::DiscountingBondEngine, ts::TermStructure) = DiscountingBondEngine(ts)
clone(pe::DiscountingSwapEngine, ts::TermStructure) = DiscountingSwapEngine(ts)
