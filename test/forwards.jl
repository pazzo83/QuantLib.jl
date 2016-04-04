using Base.Test
using QuantLib

euribor3m = euribor_index(QuantLib.Time.TenorPeriod(Dates.Month(3)))

todaysDate = Date(2006, 5, 23)
set_eval_date!(settings, todaysDate)
calendar = euribor3m.fixingCalendar
fixingDays = euribor3m.fixingDays
settlementDate = QuantLib.Time.advance(Dates.Day(fixingDays), calendar, todaysDate)

threeMonthFraQuote = zeros(10)
threeMonthFraQuote[2] = 0.030
threeMonthFraQuote[3] = 0.031
threeMonthFraQuote[4] = 0.032
threeMonthFraQuote[7] = 0.033
threeMonthFraQuote[10] = 0.034

## QUOTES ##
h1x4 = Quote(threeMonthFraQuote[2])
h2x5 = Quote(threeMonthFraQuote[3])
h3x6 = Quote(threeMonthFraQuote[4])
h6x9 = Quote(threeMonthFraQuote[7])
h9x12 = Quote(threeMonthFraQuote[10])

## RATE HELPERS ##
fraDayCount = euribor3m.dc
convention = euribor3m.convention
endOfMonth = euribor3m.endOfMonth

fra1x4 = FraRateHelper(h1x4, 1, 4, fixingDays, calendar, convention, endOfMonth, fraDayCount)
fra2x5 = FraRateHelper(h2x5, 2, 5, fixingDays, calendar, convention, endOfMonth, fraDayCount)
fra3x6 = FraRateHelper(h3x6, 3, 6, fixingDays, calendar, convention, endOfMonth, fraDayCount)
fra6x9 = FraRateHelper(h6x9, 6, 9, fixingDays, calendar, convention, endOfMonth, fraDayCount)
fra9x12 = FraRateHelper(h9x12, 9, 12, fixingDays, calendar, convention, endOfMonth, fraDayCount)

fraInstruments = FraRateHelper[fra1x4, fra2x5, fra3x6, fra6x9, fra9x12]
termStructureDC = QuantLib.Time.ISDAActualActual()
tol = 1.0e-15

fraTermStructure = PiecewiseYieldCurve(settlementDate, fraInstruments, termStructureDC, QuantLib.Math.LogLinear(), Discount(), tol, IterativeBootstrap())

euribor3m = update_termstructure(euribor3m, fraTermStructure)

fraCalendar = euribor3m.fixingCalendar
fraBusinessDayConvention = euribor3m.convention
fraFwdType = LongPosition()
fraNotional = 100.0
monthsToStart = [1, 2, 3, 6, 9]
fraTermMonths = 3
fraValueDate = QuantLib.Time.advance(Dates.Month(monthsToStart[1]), fraCalendar, settlementDate, fraBusinessDayConvention)
fraMaturityDate = QuantLib.Time.advance(Dates.Month(fraTermMonths), fraCalendar, fraValueDate, fraBusinessDayConvention)
fraStrikeRate = threeMonthFraQuote[monthsToStart[1] + 1]
myFRA = ForwardRateAgreement(fraValueDate, fraMaturityDate, fraFwdType, fraStrikeRate, fraNotional, euribor3m, fraTermStructure)

# test calculations
@test forward_rate(myFRA).rate == 0.030000000000000172
@test spot_value(myFRA) == 99.73470289713767
@test forward_value(myFRA) == 100.76666666666667
@test implied_yield(myFRA, spot_value(myFRA), forward_value(myFRA), settlementDate, SimpleCompounding(), fraDayCount).rate == 0.030039933543589872

# other tests
@test maturity_date(fra1x4) == Date(2006, 9, 26)
@test QuantLib.implied_quote(fraTermStructure.instruments[1]) == 0.030000000000000172
