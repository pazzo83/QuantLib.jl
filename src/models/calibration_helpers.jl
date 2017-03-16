type RelativePriceError <: CalibrationErrorType end
type PriceError <: CalibrationErrorType end
type ImpliedVolError <: CalibrationErrorType end

type NaiveBasketType <: CalibrationBasketType end
type MaturityStrikeByDeltaGammaBasketType <: CalibrationBasketType end

type CalibrationHelperCommon{CE <: CalibrationErrorType}
  marketValue::Float64
  calibrationErrorType::CE

  # CalibrationHelperCommon() = new{RelativePriceError}(0.0, RelativePriceError())
end

CalibrationHelperCommon() = CalibrationHelperCommon{RelativePriceError}(0.0, RelativePriceError())

type ImpliedVolatilityHelper{C <: CalibrationHelper} <: Function
  helper::C
  value::Float64
end

(iv::ImpliedVolatilityHelper)(x::Float64) = iv.value - black_price!(iv.helper, x)

function operator(iv::ImpliedVolatilityHelper)
  function _inner(x::Float64)
    return iv.value - black_price!(iv.helper, x)
  end

  return _inner
end

calibration_error{C <: CalibrationHelper}(::RelativePriceError, helper::C) =
  abs(helper.calibCommon.marketValue - model_value!(helper)) / helper.calibCommon.marketValue

# function calibration_error{C <: CalibrationHelper}(::RelativePriceError, helper::C)
#   @code_warntype model_value!(helper)
#   error("DIE")
#   return abs(helper.calibCommon.marketValue - model_value!(helper)) / helper.calibCommon.marketValue
# end

function implied_volatility!(ch::CalibrationHelper, targetValue::Float64, accuracy::Float64, maxEvals::Int, minVol::Float64, maxVol::Float64)
  ivh = ImpliedVolatilityHelper(ch, targetValue)

  solv = BrentSolver(maxEvals)

  return solve(solv, ivh, accuracy, ch.volatility.value, minVol, maxVol)
end

function update_pricing_engine{C <: CalibrationHelper, P <: PricingEngine}(ch::C, pe::P)
  newCh = clone(ch, pe)

  return newCh
end

type SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, CE <: CalibrationErrorType} <: CalibrationHelper
  lazyMixin::LazyMixin
  exerciseDate::Date
  endDate::Date
  maturity::Dm
  swapLength::Dl
  volatility::Quote
  iborIndex::IborIndex
  fixedLegTenor::TP
  fixedLegDayCount::DC_fix
  floatingLegDayCount::DC_float
  strike::Float64
  nominal::Float64
  shift::Float64
  exerciseRate::Float64
  calibCommon::CalibrationHelperCommon{CE}
  yts::T
  pricingEngine::P
  swaption::Swaption{EuropeanExercise, SettlementPhysical, P}

  SwaptionHelper(lazyMixin::LazyMixin, exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TP,
                fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64,
                calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P) =
                new{Dm, Dl, TP, DC_fix, DC_float, T, P, CE}(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor,
                fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe)

  SwaptionHelper(lazyMixin::LazyMixin, exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TP,
                fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64,
                calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P, swaption::Swaption) =
                new{Dm, Dl, TP, DC_fix, DC_float, T, P, CE}(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor,
                fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe, swaption)
end

SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine}(maturity::Dm,
              swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TP, fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float,
              yts::T, pe::P, strike::Float64 = -1.0, nominal::Float64 = 1.0, shift::Float64 = 0.0, exerciseRate::Float64 = 0.0) =
              SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, RelativePriceError}(LazyMixin(), Date(), Date(), maturity, swapLength, volatility, iborIndex, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike,
              nominal, shift, exerciseRate, CalibrationHelperCommon(), yts, pe)

SwaptionHelper{TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine}(expiryDate::Date,
              endDate::Date, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TP, fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float,
              yts::T, pe::P = NullSwaptionEngine(), strike::Float64 = -1.0, nominal::Float64 = 1.0, shift::Float64 = 0.0, exerciseRate::Float64 = 0.0) =
              SwaptionHelper{Dates.Day, Dates.Day, TP, DC_fix, DC_float, T, P, RelativePriceError}(LazyMixin(), expiryDate, endDate, Dates.Day(0), Dates.Day(0), volatility, iborIndex, fixedLegTenor,
              fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, CalibrationHelperCommon(), yts, pe)

SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, CE <: CalibrationErrorType}(lazyMixin::LazyMixin,
              exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TP,
              fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64,
              calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P, swaption::Swaption) =
              SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE}(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor,
              fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe, swaption)

SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, CE <: CalibrationErrorType}(lazyMixin::LazyMixin,
              exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex, fixedLegTenor::TP,
              fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64,
              calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P) =
              SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE}(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor,
              fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe)

function perform_calculations!(swaptionHelper::SwaptionHelper)
  calendar = swaptionHelper.iborIndex.fixingCalendar
  fixingDays = swaptionHelper.iborIndex.fixingDays
  convention = swaptionHelper.iborIndex.convention

  exerciseDate = swaptionHelper.exerciseDate == Date() ? advance(swaptionHelper.maturity, calendar, swaptionHelper.yts.referenceDate, convention) : swaptionHelper.exerciseDate

  startDate = advance(Dates.Day(fixingDays), calendar, exerciseDate, convention)

  endDate = swaptionHelper.endDate == Date() ? advance(swaptionHelper.swapLength, calendar, startDate, convention) : swaptionHelper.endDate

  fixedSchedule = QuantLib.Time.Schedule(startDate, endDate, swaptionHelper.fixedLegTenor, convention, convention, QuantLib.Time.DateGenerationForwards(), false, calendar)
  floatSchedule = QuantLib.Time.Schedule(startDate, endDate, swaptionHelper.iborIndex.tenor, convention, convention, QuantLib.Time.DateGenerationForwards(), false, calendar)

  swapEngine = DiscountingSwapEngine(swaptionHelper.yts, false)

  swapT = Receiver()

  tempSwap = VanillaSwap(swapT, swaptionHelper.nominal, fixedSchedule, 0.0, swaptionHelper.fixedLegDayCount, swaptionHelper.iborIndex, 0.0, floatSchedule, swaptionHelper.floatingLegDayCount, swapEngine)

  forward = fair_rate(tempSwap)

  if swaptionHelper.strike == -1.0
    swaptionHelper.exerciseRate = forward
  else
    swaptionHelper.exerciseRate = strike
    swapT = strike <= forward ? Receiver() : Payer()
  end

  swap = VanillaSwap(swapT, swaptionHelper.nominal, fixedSchedule, swaptionHelper.exerciseRate, swaptionHelper.fixedLegDayCount, swaptionHelper.iborIndex, 0.0, floatSchedule, swaptionHelper.floatingLegDayCount, swapEngine)
  exercise = EuropeanExercise(exerciseDate)

  swaptionHelper.swaption = Swaption(swap, exercise, swaptionHelper.pricingEngine)

  # calibration calc
  _calibration_calculate!(swaptionHelper)

  return swaptionHelper
end

function _calibration_calculate!(swaptionHelper::SwaptionHelper)
  swaptionHelper.calibCommon.marketValue = black_price!(swaptionHelper, swaptionHelper.volatility.value)

  return swaptionHelper
end

function add_times_to!(swaptionHelper::SwaptionHelper, times::Vector{Float64})
  calculate!(swaptionHelper)
  discretizedSwap = DiscretizedSwaption(swaptionHelper.swaption, reference_date(swaptionHelper.yts), swaptionHelper.yts.dc, NullLattice())
  times = vcat(times, mandatory_times(discretizedSwap))
  return times
end

function model_value!(sh::SwaptionHelper)
  calculate!(sh)
  # sh.swaption = update_pricing_engine(sh.swaption, sh.pricingEngine) # this might clone swaption
  # need to force recalc
  calculated!(sh.swaption, false)
  return npv(sh.swaption)
end

function black_price!(swaptionHelper::SwaptionHelper, sigma::Float64)
  calculate!(swaptionHelper)
  # stuff
  black = BlackSwaptionEngine(swaptionHelper.yts, Quote(sigma), Actual365(), swaptionHelper.shift)
  tempSwaption = update_pricing_engine(swaptionHelper.swaption, black)
  # swaptionHelper.swaption.pricingEngine = black
  value = npv(tempSwaption)
  return value
end

function underlying_swap!(swaptionHelper::SwaptionHelper)
  calculate!(swaptionHelper)
  return swaptionHelper.swaption.swap
end

# clone functions #
function clone(swaptionHelper::SwaptionHelper, pe::PricingEngine = swaptionHelper.pricingEngine)
  lazyMixin, calibCommon = pe == swaptionHelper.pricingEngine ? (swaptionHelper.lazyMixin, swaptionHelper.calibCommon) : (LazyMixin(), CalibrationHelperCommon())
  if (isdefined(swaptionHelper, :swaption) && swaptionHelper.pricingEngine == pe)
    newSwaptionHelper = SwaptionHelper(lazyMixin, swaptionHelper.exerciseDate, swaptionHelper.endDate, swaptionHelper.maturity, swaptionHelper.swapLength,
                        swaptionHelper.volatility, swaptionHelper.iborIndex, swaptionHelper.fixedLegTenor, swaptionHelper.fixedLegDayCount, swaptionHelper.floatingLegDayCount,
                        swaptionHelper.strike, swaptionHelper.nominal, swaptionHelper.shift, swaptionHelper.exerciseRate, calibCommon, swaptionHelper.yts,
                        pe, swaptionHelper.swaption)
  else
    newSwaptionHelper = SwaptionHelper(lazyMixin, swaptionHelper.exerciseDate, swaptionHelper.endDate, swaptionHelper.maturity, swaptionHelper.swapLength,
                        swaptionHelper.volatility, swaptionHelper.iborIndex, swaptionHelper.fixedLegTenor, swaptionHelper.fixedLegDayCount, swaptionHelper.floatingLegDayCount,
                        swaptionHelper.strike, swaptionHelper.nominal, swaptionHelper.shift, swaptionHelper.exerciseRate, calibCommon, swaptionHelper.yts,
                        pe)
  end

  return newSwaptionHelper
end
