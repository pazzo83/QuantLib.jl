using Base.Test
using QuantLib

# build a fixed rate bond to test cash flows
settlement_date = Date(2008, 9, 18)
set_eval_date!(settings, settlement_date - Dates.Day(3))
bbCurveRate = 0.055
bbDayCount = QuantLib.Time.ActualActualBond()
bbIR = InterestRate(bbCurveRate, bbDayCount, CompoundedCompounding(), QuantLib.Time.Semiannual())
ts = FlatForwardTermStructure(settlement_date, bbIR.rate, bbIR.dc, bbIR.comp, bbIR.freq)
settlement_days = 3
face_amount = 100.0
fixed_schedule = QuantLib.Time.Schedule(Date(2007, 5, 15), Date(2017, 5, 15),
                       QuantLib.Time.TenorPeriod(QuantLib.Time.Semiannual()), QuantLib.Time.Unadjusted(),
                       QuantLib.Time.Unadjusted(), QuantLib.Time.DateGenerationBackwards(), false,
                       QuantLib.Time.USGovernmentBondCalendar())

pe = DiscountingBondEngine(ts)
fixedrate_bond = FixedRateBond(settlement_days, face_amount, fixed_schedule, 0.045,
                        QuantLib.Time.ISMAActualActual(), QuantLib.Time.ModifiedFollowing(), 100.0,
                        Date(2007, 5, 15), fixed_schedule.cal, pe)

@test notional(fixedrate_bond, Date(2016, 4, 3)) == 100.0
@test accrued_amount(fixedrate_bond, Date(2016, 4, 3)) == 1.7307692307692246
@test maturity_date(fixedrate_bond) == Date(2017, 5, 15)
@test QuantLib.yield(fixedrate_bond, bbIR.dc, bbIR.comp, bbIR.freq) == 0.05506693320721388
@test duration(fixedrate_bond, bbIR, ModifiedDuration(), bbIR.dc, QuantLib.settlement_date(fixedrate_bond)) == 6.899980120600781
@test npv(fixedrate_bond) == 94.67245111530346
@test dirty_price(fixedrate_bond) == 94.67245111530346
@test clean_price(fixedrate_bond) == 93.13169024573823

@test QuantLib.get_redemption(fixedrate_bond) == fixedrate_bond.cashflows[end]
@test QuantLib.get_frequency(fixedrate_bond) == QuantLib.Time.Semiannual()
