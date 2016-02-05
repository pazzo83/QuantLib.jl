type SwapRateHelper{PrT <: Dates.Period, PrS <: Dates.Period} <: RateHelper
  rate::Quote
  tenor::PrT
  fwdStart::PrS
  swap::VanillaSwap
end

function SwapRateHelper{PrT <: Dates.Period, C <: BusinessCalendar, F <: Frequency, B <: BusinessDayConvention, DC <: DayCount, PrS <: Dates.Period, P <: PricingEngine, I <: Integer, ST <: SwapType}(rate::Float64, tenor::PrT, cal::C,
                    fixedFrequency::F, fixedConvention::B, fixedDayCount::DC, iborIndex::IborIndex, spread::Float64, fwdStart::PrS,
                    pricingEngine::P = DiscountingSwapEngine(), settlementDays::I = iborIndex.fixingDays, nominal::Float64 = 1.0, swapT::ST = Payer(), fixedRate::Float64 = 0.0)
  # do stuff
  fixedCal = cal
  floatingCal = cal
  floatTenor = iborIndex.tenor
  fixedTenor = JQuantLib.Time.TenorPeriod(fixedFrequency)
  fixedTermConvention = fixedConvention
  floatConvention = iborIndex.convention
  floatTermConvention = iborIndex.convention
  fixedRule = DateGenerationBackwards()
  floatRule = DateGenerationBackwards()
  floatDayCount = iborIndex.dc
  # fixed_rate = 0.0

  ref_date = adjust(floatingCal, floatConvention, settings.evaluation_date)
  spot_date = advance(Base.Dates.Day(settlementDays), floatingCal, ref_date, floatConvention)
  start_date = adjust(floatingCal, floatConvention, spot_date + fwdStart)
  ## TODO Float end of month (defaults to false)
  end_date = start_date + tenor

  # build schedules
  fixed_schedule = Schedule(start_date, end_date, fixedTenor, fixedConvention, fixedTermConvention, fixedRule, false, fixedCal)
  float_schedule = Schedule(start_date, end_date, floatTenor, floatConvention, floatTermConvention, floatRule, false, floatingCal)


  swap = VanillaSwap(swapT, nominal, fixed_schedule, fixedRate, fixedDayCount, iborIndex, spread, float_schedule, floatDayCount, pricingEngine, fixedConvention)

  return SwapRateHelper(Quote(rate), tenor, fwdStart, swap)
end

maturity_date(sh::SwapRateHelper) = maturity_date(sh.swap)

type DepositRateHelper <: RateHelper
  rate::Quote
  tenor::TenorPeriod
  fixingDays::Integer
  calendar::BusinessCalendar
  convention::BusinessDayConvention
  endOfMonth::Bool
  dc::DayCount
  iborIndex::IborIndex
  evaluationDate::Date
  referenceDate::Date
  earliestDate::Date
  maturityDate::Date
  fixingDate::Date
end

function DepositRateHelper(rate::Quote, tenor::TenorPeriod, fixingDays::Integer, calendar::BusinessCalendar, convention::BusinessDayConvention, endOfMonth::Bool, dc::DayCount)
  ibor_index = IborIndex("no-fix", tenor, fixingDays, NullCurrency(), calendar, convention, endOfMonth, dc)
  evaluation_date = settings.evaluation_date
  reference_date = adjust(calendar, convention, evaluation_date)
  earliest_date = advance(Dates.Day(fixingDays), calendar, reference_date, convention)
  maturity_d = maturity_date(ibor_index, earliest_date)
  fix_date = fixing_date(ibor_index, earliest_date)
  return DepositRateHelper(rate, tenor, fixingDays, calendar, convention, endOfMonth, dc, ibor_index, evaluation_date, reference_date, earliest_date, maturity_d, fix_date)
end

maturity_date{R <: RateHelper}(rate::R) = rate.maturityDate

value{R <: RateHelper}(rate::R) = rate.rate.value

function implied_quote(depo::DepositRateHelper)
  return fixing(depo.iborIndex, depo.iborIndex.ts, depo.fixingDate, true)
end

function implied_quote(swap_helper::SwapRateHelper)
  swap = swap_helper.swap
  recalculate!(swap)
  #
  # println("Floating Leg NPV ", floating_leg_NPV(swap))
  # println("Floating Leg BPS ", floating_leg_BPS(swap))
  # println("Fixed Leg BPS ", fixed_leg_BPS(swap))
  # println("Swap spread ", swap.spread)

  floatingLegNPV = floating_leg_NPV(swap)
  spread = swap.spread
  spreadNPV = floating_leg_BPS(swap) / basisPoint * spread
  totNPV = -(floatingLegNPV + spreadNPV)

  return totNPV / (fixed_leg_BPS(swap) / basisPoint)
end
