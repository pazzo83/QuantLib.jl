type MidPointCdsEngine{P <: AbstractDefaultProbabilityTermStructure, Y <: YieldTermStructure} <: PricingEngine
  probability::P
  recoveryRate::Float64
  discountCurve::Y
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
