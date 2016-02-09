type SpreadCDSHelper{P <: Dates.Period, I <: Integer, C <: BusinessCalendar, F <: Frequency, PC <: BusinessDayConvention, DC <: DayCount, YTS <: YieldTermStructure} <: AbstractCDSHelper
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
  # swap::CreditDefaultSwap

  # main constructor
  # SpreadCDSHelper(runningSpread::Quote, tenor::P, settlementDays::I, calendar::C, frequency::F, paymentConvention::PC, dc::DC, recoveryRate::Float64, schedule::Schedule,
  #                 discountCurve::YTS, settlesAccrual::Bool, paysAtDefaultTime::Bool, protectionStart::Date) =
  #                 new(runningSpread, tenor, settlementDays, calendar, frequency, paymentConvention, dc, recoveryRate, schedule, discountCurve, settlesAccrual, paysAtDefaultTime, protectionStart)
end

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
