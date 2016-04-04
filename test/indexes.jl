using Base.Dates
using QuantLib

ts = FlatForwardTermStructure(settlement_date, bbIR.rate, bbIR.dc, bbIR.comp, bbIR.freq)

# Construct a 6 month euribor index
euribor6m = euribor_index(QuantLib.Time.TenorPeriod(Base.Dates.Month(6)))

# Construct a 3 month usd libor index
libor_3m = usd_libor_index(QuantLib.Time.TenorPeriod(Base.Dates.Month(3)), ts)

# test date calculations
fixing_date(euribor6m, Date(2016, 4, 3)) == Date(2016, 3, 31)
maturity_date(euribor6m, Date(2016, 4, 3)) == Date(2016, 10, 3)
QuantLib.value_date(euribor6m, Date(2016, 4, 3)) == Date(2016, 4, 5)

fixing_date(libor_3m, Date(2016, 4, 3)) == Date(2016, 3, 31)
maturity_date(libor_3m, Date(2016, 4, 3)) == Date(2016, 7, 5)
QuantLib.value_date(libor_3m, Date(2016, 4, 3)) == Date(2016, 4, 5)

# other calculations
fixing(euribor6m, ts, fixing_date(euribor6m, Date(2016, 4, 3))) == 0.0540983606557383
fixing(libor_3m, ts, fixing_date(libor_3m, Date(2016, 4, 3))) == 0.05343944335694899
