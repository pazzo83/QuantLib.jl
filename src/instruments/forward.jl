type ForwardRateAgreement{DC <: DayCount, BC <: BusinessCalendar, C <: BusinessDayConvention, Y <: YieldTermStructure, P <: PositionType} <: AbstractForward
  lazyMixin::LazyMixin
  underlyingIncome::Float64
  underlyingSpotValue::Float64
  dc::DC
  calendar::BC
  convention::C
  settlementDays::Int
  payoff::ForwardTypePayoff
  valueDate::Date
  maturityDate::Date
  discountCurve::Y
  fraType::P
  forwardRate::InterestRate
  strikeForwardRate::InterestRate
  notionalAmount::Float64
  iborIndex::IborIndex
end

function ForwardRateAgreement(valueDate::Date, maturityDate::Date, position::PositionType, strikeForward::Float64, notionalAmount::Float64,
                              iborIndex::IborIndex, discountCurve::YieldTermStructure)
  calendar = iborIndex.fixingCalendar
  convention = iborIndex.convention
  settlementDays = iborIndex.fixingDays
  maturityDate = adjust(calendar, maturityDate, convention)
  fixingDate = advance(Dates.Day(-settlementDays), calendar, valueDate)
  forwardRate = InterestRate(fixing(iborIndex, iborIndex.ts, fixingDate), iborIndex.dc, SimpleCompounding(), Once())
  strikeForwardRate = InterestRate(strikeForward, iborIndex.dc, SimpleCompounding(), Once())

  strike = notionalAmount * compound_factor(strikeForwardRate, valueDate, maturityDate)
  payoff = ForwardTypePayoff(position, strike)

  return ForwardRateAgreement(LazyMixin(), 0.0, 0.0, iborIndex.dc, calendar, convention, settlementDays, payoff, valueDate, maturityDate, discountCurve,
                              position, forwardRate, strikeForwardRate, notionalAmount, iborIndex)
end
