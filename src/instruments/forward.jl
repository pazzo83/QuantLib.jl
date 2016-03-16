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
  maturityDate = adjust(calendar, convention, maturityDate)
  fixingDate = advance(Dates.Day(-settlementDays), calendar, valueDate)
  forwardRate = InterestRate(fixing(iborIndex, iborIndex.ts, fixingDate), iborIndex.dc, SimpleCompounding(), Once())
  strikeForwardRate = InterestRate(strikeForward, iborIndex.dc, SimpleCompounding(), Once())

  strike = notionalAmount * compound_factor(strikeForwardRate, valueDate, maturityDate)
  payoff = ForwardTypePayoff(position, strike)

  return ForwardRateAgreement(LazyMixin(), 0.0, 0.0, iborIndex.dc, calendar, convention, settlementDays, payoff, valueDate, maturityDate, discountCurve,
                              position, forwardRate, strikeForwardRate, notionalAmount, iborIndex)
end

function spot_value(fra::ForwardRateAgreement)
  calculate!(fra)
  return fra.notionalAmount * compound_factor(forward_rate(fra), fra.valueDate, fra.maturityDate) * discount(fra.discountCurve, fra.maturityDate)
end

function forward_rate(fra::ForwardRateAgreement)
  calculate!(fra)
  return fra.forwardRate
end

function forward_value(fra::ForwardRateAgreement)
  calculate!(fra)
  return (fra.underlyingSpotValue - fra.underlyingIncome) / discount(fra.discountCurve, fra.maturityDate)
end

spot_income(::ForwardRateAgreement, ::YieldTermStructure) = 0.0

function implied_yield(fra::ForwardRateAgreement, underlyingSpotValue::Float64, forwardValue::Float64, settlementDate::Date, comp::CompoundingType, dc::DayCount)
  t = year_fraction(dc, settlementDate, fra.maturityDate)
  compoundingFactor = forwardValue / (underlyingSpotValue - spot_income(fra, fra.discountCurve))

  return implied_rate(compoundingFactor, dc, comp, t, Annual())
end

function perform_calculations!(fra::ForwardRateAgreement)
  fixingDate = advance(Dates.Day(-fra.settlementDays), fra.calendar, fra.valueDate)
  fra.forwardRate = InterestRate(fixing(fra.iborIndex, fra.iborIndex.ts, fixingDate), fra.iborIndex.dc, SimpleCompounding(), Once())
  fra.underlyingSpotValue = spot_value(fra)
  fra.underlyingIncome = 0.0
  # NPV calc

  return fra
end
