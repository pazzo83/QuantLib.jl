immutable SwapIndex{S <: AbstractString, TP <: TenorPeriod, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure} <: InterestRateIndex
  familyName::S
  tenor::TP
  fixingDays::Int
  currency::Currency
  fixingCalendar::B
  fixedLegTenor::Dates.Period
  fixedLegConvention::C
  fixedLegDayCount::DC
  discount::T
  iborIndex::IborIndex
  exogenousDiscount::Bool
  lastSwap::VanillaSwap
  lastFixingDate::Date

  SwapIndex(familyName::S, tenor::TP, fixingDays::Int, currency::Currency, fixingCalendar::B, fixedLegTenor::Dates.Period, fixedLegConvention::C, fixedLegDayCount::DC,
            discount::T, iborIndex::IborIndex, exogenousDiscount::Bool = true) =
            new{S, TP, B, C, DC, T}(familyName, tenor, fixingDays, currency, fixingCalendar, fixedLegTenor, fixedLegConvention, fixedLegDayCount, discount, iborIndex, exogenousDiscount)
end

SwapIndex{S <: AbstractString, TP <: TenorPeriod, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure}(familyName::S, tenor::TP, settlementDays::Int, currency::Currency, fixingCalendar::B,
          fixedLegTenor::Dates.Period, fixedLegConvention::C, fixedLegDayCount::DC, discount::T, iborIndex::IborIndex) = SwapIndex{S, TP, B, C, DC, T}(familyName, tenor, settlementDays, currency, fixingCalendar, fixedLegTenor,
          fixedLegConvention, fixedLegDayCount, discount, iborIndex)

function EuriborSwapIsdaFixA(tenor::TenorPeriod, forwardingTS::YieldTermStructure, discTS::YieldTermStructure)
  ibor = Int(tenor.period) > 1 ? euribor_index(TenorPeriod(Dates.Month(6)), forwardingTS) : euribor_index(TenorPeriod(Dates.Month(3)), forwardingTS)

  return SwapIndex("EuriborSwapIsdaFixA", tenor, 2, EURCurrency(), TargetCalendar(), Dates.Year(1), ModifiedFollowing(), BondThirty360(), discTS, ibor)
end
