struct BlackSwaptionEngine{Y <: YieldTermStructure, S <: SwaptionVolatilityStructure, DC <: DayCount} <: PricingEngine
  yts::Y
  vol::Quote
  volStructure::S
  dc::DC
  displacement::Float64
end

function BlackSwaptionEngine(yts::Y, vol::Quote, dc::DC, displacement::Float64 = 0.0) where {Y <: YieldTermStructure, DC <: DayCount}
  cwv = ConstantSwaptionVolatility(0, QuantLib.Time.NullCalendar(), QuantLib.Time.Following(), vol, dc)
  BlackSwaptionEngine{Y, typeof(cwv), DC}(yts, vol, cwv, dc, displacement)
end

# pricing methods #
function black_formula(optionType::OptionType, strike::Float64, forward::Float64, stdDev::Float64, discount::Float64 = 1.0, displacement::Float64 = 0.0)
  # TODO check requirements (see cpp)
  opt_type = QuantLib.value(optionType)
  if stdDev == 0.0
    return max((forward - strike) * opt_type, 0.0 * discount)
  end

  forward += displacement
  strike += displacement

  if strike == 0.0
    return isa(optionType, Call) ? forward * discount : 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev
  d2 = d1 - stdDev
  norm = Normal() # using distributions.jl
  nd1 = cdf(norm, opt_type * d1)
  nd2 = cdf(norm, opt_type * d2)
  result = discount * opt_type * (forward * nd1 - strike * nd2)

  return result
end

function black_formula_standard_dev_derivative(strike::Float64, forward::Float64, stdDev::Float64, discount::Float64, displacement::Float64)
  forward += displacement
  strike += displacement

  if stdDev == 0.0 || strike == 0.0
    return 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev

  return discount * forward * distribution_derivative(Normal(), d1)
end

get_annuity(delivery::SettlementPhysical, swap::VanillaSwap) = abs(fixed_leg_BPS(swap)) / basisPoint

# calculation methods #
function _calculate!(pe::BlackSwaptionEngine, swaption::Swaption)
  exerciseDate = swaption.exercise.dates[1]
  swap = swaption.swap

  strike = swap.fixedRate

  # override swap's pricing engine temporarily, bypassing normal calc flow, since swap.iborIndex might be using a diff curve
  tempSwap = clone(swap, DiscountingSwapEngine{typeof(pe.yts)}(pe.yts))
  # _calculate!(DiscountingSwapEngine(pe.yts), swap)
  calculate!(tempSwap)
  # swap.lazyMixin.calculated = true
  atmForward = fair_rate(tempSwap)

  if tempSwap.spread != 0.0
    correction = tempSwap.spread * abs(floating_leg_BPS(tempSwap) / fixed_leg_BPS(tempSwap))
    strike -= correction
    atmForward -= correction
    swaption.results.additionalResults["spreadCorrection"] = correction
  else
    swaption.results.additionalResults["spreadCorrection"] = 0.0
  end

  swaption.results.additionalResults["strike"] = strike
  swaption.results.additionalResults["atmForward"] = atmForward

  annuity = get_annuity(swaption.delivery, tempSwap)
  swaption.results.additionalResults["annuity"] = annuity

  # the swap length calculation might be improved using the value date of the exercise date
  swapLength = swap_length(pe.volStructure, exerciseDate, date(get_latest_coupon(tempSwap.legs[1])))
  swaption.results.additionalResults["swapLength"] = swapLength

  variance = black_variance(pe.volStructure, exerciseDate, swapLength, strike)

  stdDev = sqrt(variance)
  swaption.results.additionalResults["stdDev"] = stdDev
  w = isa(tempSwap.swapT, Payer) ? Call() : Put()

  swaption.results.value = black_formula(w, strike, atmForward, stdDev, annuity, pe.displacement)

  exerciseTime = time_from_reference(pe.volStructure, exerciseDate)

  swaption.results.additionalResults["vega"] = sqrt(exerciseTime) * black_formula_standard_dev_derivative(strike, atmForward, stdDev, annuity, pe.displacement)

  # # resetting swap
  # reset!(swap.results)
  # swap.lazyMixin.calculated = false

  return swaption
end
