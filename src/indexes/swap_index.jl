struct SwapIndex{TP <: TenorPeriod, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure, II <: IborIndex} <: InterestRateIndex
  familyName::String
  tenor::TP
  fixingDays::Int
  currency::Currency
  fixingCalendar::B
  fixedLegTenor::Dates.Period
  fixedLegConvention::C
  fixedLegDayCount::DC
  discount::T
  iborIndex::II
  exogenousDiscount::Bool
  lastSwap::VanillaSwap
  lastFixingDate::Date

  SwapIndex{TP, B, C, DC, T, II}(familyName::String,
                                tenor::TP,
                                fixingDays::Int,
                                currency::Currency,
                                fixingCalendar::B,
                                fixedLegTenor::Dates.Period,
                                fixedLegConvention::C,
                                fixedLegDayCount::DC,
                                discount::T,
                                iborIndex::II,
                                exogenousDiscount::Bool = true) where {TP, B, C, DC, T, II} =
            new{TP, B, C, DC, T, II}(familyName, tenor, fixingDays, currency, fixingCalendar, fixedLegTenor, fixedLegConvention, fixedLegDayCount, discount,
                                    iborIndex, exogenousDiscount)
end

SwapIndex(familyName::String,
          tenor::TP,
          settlementDays::Int,
          currency::Currency,
          fixingCalendar::B,
          fixedLegTenor::Dates.Period,
          fixedLegConvention::C,
          fixedLegDayCount::DC, discount::T,
          iborIndex::II) where {TP <: TenorPeriod, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure, II <: IborIndex} =
          SwapIndex{TP, B, C, DC, T, }(familyName, tenor, settlementDays, currency, fixingCalendar, fixedLegTenor, fixedLegConvention, fixedLegDayCount, discount, iborIndex)

function EuriborSwapIsdaFixA(tenor::TenorPeriod, forwardingTS::YieldTermStructure, discTS::YieldTermStructure)
  ibor = tenor.period.value > 1 ? euribor_index(TenorPeriod(Dates.Month(6)), forwardingTS) : euribor_index(TenorPeriod(Dates.Month(3)), forwardingTS)

  return SwapIndex("EuriborSwapIsdaFixA", tenor, 2, EURCurrency(), TargetCalendar(), Dates.Year(1), ModifiedFollowing(), BondThirty360(), discTS, ibor)
end
