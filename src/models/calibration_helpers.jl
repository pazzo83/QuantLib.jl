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

type SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, CE <: CalibrationErrorType, ST <: SwapType, B <: BusinessDayConvention, SP <: PricingEngine, STP <: TenorPeriod, CUR <: AbstractCurrency, IB <: BusinessCalendar, IC <: BusinessDayConvention, IDC <: DayCount, IT <: TermStructure} <: CalibrationHelper
  lazyMixin::LazyMixin
  exerciseDate::Date
  endDate::Date
  maturity::Dm
  swapLength::Dl
  volatility::Quote
  iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}
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
  swaption::Swaption{EuropeanExercise, SettlementPhysical, P, ST, DC_fix, DC_float, B, Leg, SP, STP, CUR, IB, IC, IDC, IT}



  # SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, ST, B, SP, STP, CUR, IB, IC, IDC, IT}(lazyMixin::LazyMixin, exerciseDate::Date, endDate::Date, maturity::Dm,
  #               swapLength::Dl, volatility::Quote, iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}, fixedLegTenor::TP, fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float,
  #               strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64, calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P,
  #               swaption::Swaption{EuropeanExercise, SettlementPhysical, P, ST, DC_fix, DC_float, B, Leg, SP, STP, CUR, IB, IC, IDC, IT}) =
  #               new{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, ST, B, SP, STP, CUR, IB, IC, IDC, IT}(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex,
  #               fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe, swaption)
end

function SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, STP, CUR, IB, IC, IDC, IT}(lazyMixin::LazyMixin, exerciseDate::Date, endDate::Date, maturity::Dm,
              swapLength::Dl, volatility::Quote, iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}, fixedLegTenor::TP, fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float,
              strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64, calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P)

  swaption = build_swaption(exerciseDate, endDate, maturity, swapLength, iborIndex, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, exerciseRate, yts, pe)

  SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, typeof(swaption.swap.swapT), typeof(swaption.swap.paymentConvention), typeof(swaption.swap.pricingEngine), STP, CUR, IB, IC, IDC, IT}(lazyMixin, exerciseDate,
  endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor, fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe,
  swaption)
end

SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, STP <: TenorPeriod, CUR <: AbstractCurrency, IB <: BusinessCalendar, IC <: BusinessDayConvention, IDC <: DayCount, IT <: TermStructure}(maturity::Dm,
              swapLength::Dl, volatility::Quote, iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}, fixedLegTenor::TP, fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float,
              yts::T, pe::P, strike::Float64 = -1.0, nominal::Float64 = 1.0, shift::Float64 = 0.0, exerciseRate::Float64 = 0.0) =
              SwaptionHelper(LazyMixin(), Date(), Date(), maturity, swapLength, volatility, iborIndex, fixedLegTenor,
              fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, CalibrationHelperCommon(), yts, pe)

SwaptionHelper{TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, STP <: TenorPeriod, CUR <: AbstractCurrency, IB <: BusinessCalendar, IC <: BusinessDayConvention, IDC <: DayCount, IT <: TermStructure}(expiryDate::Date,
              endDate::Date, volatility::Quote, iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}, fixedLegTenor::TP, fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float,
              yts::T, pe::P = NullSwaptionEngine(), strike::Float64 = -1.0, nominal::Float64 = 1.0, shift::Float64 = 0.0, exerciseRate::Float64 = 0.0) =
              SwaptionHelper(LazyMixin(), expiryDate, endDate, Dates.Day(0), Dates.Day(0), volatility, iborIndex, fixedLegTenor,
              fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, CalibrationHelperCommon(), yts, pe)

# SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, CE <: CalibrationErrorType, ST <: SwapType, B <: BusinessDayConvention, SP <: PricingEngine, STP <: TenorPeriod, CUR <: AbstractCurrency, IB <: BusinessCalendar, IC <: BusinessDayConvention, IDC <: DayCount, IT <: TermStructure}(lazyMixin::LazyMixin,
#               exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}, fixedLegTenor::TP,
#               fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64,
#               calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P,
#               swaption::Swaption{EuropeanExercise, SettlementPhysical, P, ST, DC_fix, DC_float, B, Leg, SP, STP, CUR, IB, IC, IDC, IT}) =
#               SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, ST, B, SP, STP, CUR, IB, IC, IDC, IT}(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor,
#               fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe, swaption)

# SwaptionHelper{Dm <: Dates.Period, Dl <: Dates.Period, TP <: TenorPeriod, DC_fix <: DayCount, DC_float <: DayCount, T <: YieldTermStructure, P <: PricingEngine, CE <: CalibrationErrorType, STP <: TenorPeriod, CUR <: AbstractCurrency, IB <: BusinessCalendar, IC <: BusinessDayConvention, IDC <: DayCount, IT <: TermStructure}(lazyMixin::LazyMixin,
#               exerciseDate::Date, endDate::Date, maturity::Dm, swapLength::Dl, volatility::Quote, iborIndex::IborIndex{STP, CUR, IB, IC, IDC, IT}, fixedLegTenor::TP,
#               fixedLegDayCount::DC_fix, floatingLegDayCount::DC_float, strike::Float64, nominal::Float64, shift::Float64, exerciseRate::Float64,
#               calibCommon::CalibrationHelperCommon{CE}, yts::T, pe::P) =
#               SwaptionHelper(lazyMixin, exerciseDate, endDate, maturity, swapLength, volatility, iborIndex, fixedLegTenor,
#               fixedLegDayCount, floatingLegDayCount, strike, nominal, shift, exerciseRate, calibCommon, yts, pe)



function perform_calculations!(swaptionHelper::SwaptionHelper)
  # calendar = swaptionHelper.iborIndex.fixingCalendar
  # fixingDays = swaptionHelper.iborIndex.fixingDays
  # convention = swaptionHelper.iborIndex.convention
  #
  # exerciseDate = swaptionHelper.exerciseDate == Date() ? advance(swaptionHelper.maturity, calendar, swaptionHelper.yts.referenceDate, convention) : swaptionHelper.exerciseDate
  #
  # startDate = advance(Dates.Day(fixingDays), calendar, exerciseDate, convention)
  #
  # endDate = swaptionHelper.endDate == Date() ? advance(swaptionHelper.swapLength, calendar, startDate, convention) : swaptionHelper.endDate
  #
  # fixedSchedule = QuantLib.Time.Schedule(startDate, endDate, swaptionHelper.fixedLegTenor, convention, convention, QuantLib.Time.DateGenerationForwards(), false, calendar)
  # floatSchedule = QuantLib.Time.Schedule(startDate, endDate, swaptionHelper.iborIndex.tenor, convention, convention, QuantLib.Time.DateGenerationForwards(), false, calendar)
  #
  # swapEngine = DiscountingSwapEngine(swaptionHelper.yts, false)
  #
  # swapT = Receiver()
  #
  # tempSwap = VanillaSwap(swapT, swaptionHelper.nominal, fixedSchedule, 0.0, swaptionHelper.fixedLegDayCount, swaptionHelper.iborIndex, 0.0, floatSchedule, swaptionHelper.floatingLegDayCount, swapEngine)
  #
  # forward = fair_rate(tempSwap)
  #
  # if swaptionHelper.strike == -1.0
  #   swaptionHelper.exerciseRate = forward
  # else
  #   swaptionHelper.exerciseRate = strike
  #   swapT = strike <= forward ? Receiver() : Payer()
  # end
  #
  # swap = VanillaSwap(swapT, swaptionHelper.nominal, fixedSchedule, swaptionHelper.exerciseRate, swaptionHelper.fixedLegDayCount, swaptionHelper.iborIndex, 0.0, floatSchedule, swaptionHelper.floatingLegDayCount, swapEngine)
  # exercise = EuropeanExercise(exerciseDate)
  #
  # swaptionHelper.swaption = Swaption(swap, exercise, swaptionHelper.pricingEngine)

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

function build_swaption(exerciseDate::Date, endDate::Date, maturity::Dates.Period, swapLength::Dates.Period, iborIndex::IborIndex,
                        fixedLegTenor::TenorPeriod, fixedLegDayCount::DayCount, floatingLegDayCount::DayCount, strike::Float64, nominal::Float64,
                        exerciseRate::Float64, yts::TermStructure, pe::PricingEngine)
  calendar = iborIndex.fixingCalendar
  fixingDays = iborIndex.fixingDays
  convention = iborIndex.convention

  exerciseDate = exerciseDate == Date() ? advance(maturity, calendar, yts.referenceDate, convention) : exerciseDate

  startDate = advance(Dates.Day(fixingDays), calendar, exerciseDate, convention)

  endDate = endDate == Date() ? advance(swapLength, calendar, startDate, convention) : endDate

  fixedSchedule = QuantLib.Time.Schedule(startDate, endDate, fixedLegTenor, convention, convention, QuantLib.Time.DateGenerationForwards(), false, calendar)
  floatSchedule = QuantLib.Time.Schedule(startDate, endDate, iborIndex.tenor, convention, convention, QuantLib.Time.DateGenerationForwards(), false, calendar)

  swapEngine = DiscountingSwapEngine(yts, false)

  swapT = Receiver()

  tempSwap = VanillaSwap(swapT, nominal, fixedSchedule, 0.0, fixedLegDayCount, iborIndex, 0.0, floatSchedule, floatingLegDayCount, swapEngine)

  forward = fair_rate(tempSwap)

  if strike == -1.0
    exerciseRate = forward
  else
    exerciseRate = strike
    swapT = strike <= forward ? Receiver() : Payer()
  end

  swap = VanillaSwap(swapT, nominal, fixedSchedule, exerciseRate, fixedLegDayCount, iborIndex, 0.0, floatSchedule, floatingLegDayCount, swapEngine)
  exercise = EuropeanExercise(exerciseDate)

  swaption = Swaption(swap, exercise, pe)

  return swaption
end

function build_swaption(sh::SwaptionHelper, pe::PricingEngine)
  # for cloning
  return build_swaption(sh.exerciseDate, sh.endDate, sh.maturity, sh.swapLength, sh.iborIndex, sh.fixedLegTenor, sh.fixedLegDayCount, sh.floatingLegDayCount,
                        sh.strike, sh.nominal, sh.exerciseRate, sh.yts, pe)
end

# clone functions #
function clone{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, ST, B, SP, STP, CUR, IB, IC, IDC, IT, P2}(swaptionHelper::SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P, CE, ST, B, SP, STP, CUR, IB, IC, IDC, IT}, pe::P2 = swaptionHelper.pricingEngine)
  lazyMixin, calibCommon = pe == swaptionHelper.pricingEngine ? (swaptionHelper.lazyMixin, swaptionHelper.calibCommon) : (LazyMixin(), CalibrationHelperCommon())
  swaption = build_swaption(swaptionHelper, pe)
  newSwaptionHelper = SwaptionHelper{Dm, Dl, TP, DC_fix, DC_float, T, P2, CE, ST, B, SP, STP, CUR, IB, IC, IDC, IT}(lazyMixin, swaptionHelper.exerciseDate, swaptionHelper.endDate, swaptionHelper.maturity, swaptionHelper.swapLength,
                      swaptionHelper.volatility, swaptionHelper.iborIndex, swaptionHelper.fixedLegTenor, swaptionHelper.fixedLegDayCount, swaptionHelper.floatingLegDayCount,
                      swaptionHelper.strike, swaptionHelper.nominal, swaptionHelper.shift, swaptionHelper.exerciseRate, calibCommon, swaptionHelper.yts,
                      pe, swaption)


  return newSwaptionHelper
end
