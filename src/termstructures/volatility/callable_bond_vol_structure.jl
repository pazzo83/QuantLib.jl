type CallableBondConstantVol{C <: BusinessCalendar, DC <: DayCount, D <: Dates.Period} <: CallableBondVolatilityStructure
  settlementDays::Int
  calendar::C
  volatility::Quote
  dc::DC
  maxBondTenor::D
end

CallableBondConstantVol(settlementDays::Int, cal::BusinessCalendar, fwdYieldVol::Quote, dc::DayCount) =
                      CallableBondConstantVol(settlementDays, cal, fwdYieldVol, dc, Dates.Year(100))
