using Base.Test
using QuantLib

cal = QuantLib.Time.TargetCalendar()
settlementDate = Date(2002, 2, 19)
todays_date = Date(2002, 2, 15)
set_eval_date!(settings, todays_date)

flat_rate = Quote(0.04875825)
rhTermStructure = FlatForwardTermStructure(settlementDate, cal, flat_rate, QuantLib.Time.Actual365())

# Define swap
fixedLegFrequency = QuantLib.Time.Annual()
fixedLegConvention = QuantLib.Time.Unadjusted()
floatingLegConvention = QuantLib.Time.ModifiedFollowing()
fixedLegDayCounter = QuantLib.Time.EuroThirty360()
floatingLegFrequency = QuantLib.Time.Semiannual()

swapType = Payer()
dummyFixedRate = 0.03
indexSixMonths = euribor_index(QuantLib.Time.TenorPeriod(Dates.Month(6)), rhTermStructure)

startDate = QuantLib.Time.advance(Dates.Year(1), cal, settlementDate, floatingLegConvention)
maturity = QuantLib.Time.advance(Dates.Year(5), cal, startDate, floatingLegConvention)

fixedSchedule = QuantLib.Time.Schedule(startDate, maturity, QuantLib.Time.TenorPeriod(fixedLegFrequency), fixedLegConvention, fixedLegConvention, QuantLib.Time.DateGenerationForwards(), false, cal)
floatSchedule = QuantLib.Time.Schedule(startDate, maturity, QuantLib.Time.TenorPeriod(floatingLegFrequency), floatingLegConvention, floatingLegConvention, QuantLib.Time.DateGenerationForwards(), false, cal)

swap = VanillaSwap(swapType, 1000.0, fixedSchedule, dummyFixedRate, fixedLegDayCounter, indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))

# Calculations
@test fair_rate(swap) == 0.050000005910246725

# Other tests
@test QuantLib.get_fixed_reset_dates(swap) == Date[Date(2003, 2, 19), Date(2004, 2, 19), Date(2005, 2, 19), Date(2006, 2, 19), Date(2007, 2, 19)]
@test QuantLib.get_fixed_pay_dates(swap) == Date[Date(2004, 2, 19), Date(2005, 2, 21), Date(2006, 2, 20), Date(2007, 2, 19), Date(2008, 2, 19)]
@test QuantLib.get_floating_reset_dates(swap) == Date[Date(2003, 2, 19), Date(2003, 8, 19), Date(2004, 2, 19), Date(2004, 8, 19), Date(2005, 2, 21), Date(2005, 8, 19), Date(2006, 2, 20), Date(2006, 8, 21), Date(2007, 2, 19), Date(2007, 8, 20)]
@test QuantLib.get_floating_pay_dates(swap) == Date[Date(2003, 8, 19), Date(2004, 2, 19), Date(2004, 8, 19), Date(2005, 2, 21), Date(2005, 8, 19), Date(2006, 2, 20), Date(2006, 8, 21), Date(2007, 2, 19), Date(2007, 8, 20), Date(2008, 2, 19)]
@test QuantLib.get_floating_accrual_times(swap) == Float64[0.5027777777777778,0.5111111111111111,0.5055555555555555,0.5166666666666667,0.49722222222222223,0.5138888888888888,0.5055555555555555,0.5055555555555555,0.5055555555555555,0.5083333333333333]
@test QuantLib.get_floating_spreads(swap) == Float64[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]


## CDS ##
