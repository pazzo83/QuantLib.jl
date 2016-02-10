type SpreadCDSHelper{P <: Dates.Period, I <: Integer, C <: BusinessCalendar, F <: Frequency, PC <: BusinessDayConvention, DC <: DayCount, YTS <: YieldTermStructure, PROB <: AbstractDefaultProbabilityTermStructure} <: AbstractCDSHelper
  runningSpread::Quote
  tenor::P
  settlementDays::I
  calendar::C
  frequency::F
  paymentConvention::PC
  dc::DC
  recoveryRate::Float64
  schedule::Schedule
  discountCurve::YTS
  settlesAccrual::Bool
  paysAtDefaultTime::Bool
  protectionStart::Date
  probability::PROB
  swap::CreditDefaultSwap

  # main constructor
  SpreadCDSHelper(runningSpread::Quote, tenor::P, settlementDays::I, calendar::C, frequency::F, paymentConvention::PC, dc::DC, recoveryRate::Float64, schedule::Schedule,
                  discountCurve::YTS, settlesAccrual::Bool, paysAtDefaultTime::Bool, protectionStart::Date, probability::PROB) =
                  new(runningSpread, tenor, settlementDays, calendar, frequency, paymentConvention, dc, recoveryRate, schedule, discountCurve, settlesAccrual, paysAtDefaultTime, protectionStart, probability)
end

# Catch all constructor
SpreadCDSHelper{P <: Dates.Period, I <: Integer, C <: BusinessCalendar, F <: Frequency, PC <: BusinessDayConvention, DC <: DayCount, YTS <: YieldTermStructure, PROB <: AbstractDefaultProbabilityTermStructure}(
                runningSpread::Quote, tenor::P, settlementDays::I, calendar::C, frequency::F, paymentConvention::PC, dc::DC, recoveryRate::Float64,
                schedule::Schedule, discountCurve::YTS, settlesAccrual::Bool, paysAtDefaultTime::Bool, protectionStart::Date, probability::PROB = NullDefaultProbabilityTermStructure()) =

                SpreadCDSHelper{P, I, C, F, PC, DC, YTS, PROB}(runningSpread, tenor, settlementDays, calendar, frequency, paymentConvention, dc, recoveryRate, schedule, discountCurve,
                                settlesAccrual, paysAtDefaultTime, protectionStart, probability)

# main outer constructor
function SpreadCDSHelper(runningSpread::Quote, tenor::Dates.Period, settlementDays::Int, calendar::BusinessCalendar, frequency::Frequency, paymentConvention::BusinessDayConvention, rule::DateGenerationRule,
                        dc::DayCount, recoveryRate::Float64, discountCurve::YieldTermStructure = NullYieldTermStructure(), settlesAccrual::Bool = true, paysAtDefaultTime::Bool = true)
  evalDate = settings.evaluation_date
  protectionStart, schedule = initialize_dates(SpreadCDSHelper, evalDate, settlementDays, calendar, paymentConvention, tenor, frequency, rule)

  return SpreadCDSHelper(runningSpread, tenor, settlementDays, calendar, frequency, paymentConvention, dc, recoveryRate, schedule, discountCurve, settlesAccrual, paysAtDefaultTime, protectionStart)
end

function initialize_dates(::Type{SpreadCDSHelper}, evaluationDate::Date, settlementDays::Int, calendar::BusinessCalendar, paymentConvention::BusinessDayConvention, tenor::Dates.Period, frequency::Frequency, rule::DateGenerationRule)
  protectionStart = evaluationDate + Dates.Day(settlementDays)
  startDate = adjust(calendar, paymentConvention, protectionStart)
  endDate = evaluationDate + tenor

  schedule = Schedule(startDate, endDate, TenorPeriod(frequency), paymentConvention, Unadjusted(), rule, false, calendar)

  return protectionStart, schedule
end

maturity_date(cds::SpreadCDSHelper) = cds.schedule.dates[end]

function reset_engine!(cds::SpreadCDSHelper)
  cds.swap = CreditDefaultSwap(Buyer(), 100.0, 0.01, cds.schedule, cds.paymentConvention, cds.dc, cds.settlesAccrual, cds.paysAtDefaultTime,
              cds.protectionStart, MidPointCdsEngine(cds.probability, cds.recoveryRate, cds.discountCurve))

  return cds
end

# defaults to existing prob ts, note: this doesn't copy the swap!
clone(cds::SpreadCDSHelper, probTs::AbstractDefaultProbabilityTermStructure = cds.probability) =
      SpreadCDSHelper(cds.runningSpread, cds.tenor, cds.settlementDays, cds.calendar, cds.frequency, cds.paymentConvention, cds.dc, cds.recoveryRate,
                      cds.schedule, cds.discountCurve, cds.settlesAccrual, cds.paysAtDefaultTime, cds.protectionStart, probTs)
