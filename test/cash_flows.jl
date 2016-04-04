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

# generate a floating rate bond
floating_cal = QuantLib.Time.TargetCalendar()
settlement_days = 3
cap_vol = ConstantOptionVolatility(3, floating_cal, QuantLib.Time.ModifiedFollowing(), 0.0, QuantLib.Time.Actual365())

fb_issue_date = Date(2005, 10, 21)
fb_bond_engine = DiscountingBondEngine(ts)
fb_dc = QuantLib.Time.Actual360()
conv = QuantLib.Time.ModifiedFollowing()

fb_schedule = QuantLib.Time.Schedule(Date(2005, 10, 21), Date(2010, 10, 21), QuantLib.Time.TenorPeriod(QuantLib.Time.Quarterly()),
                                        QuantLib.Time.Unadjusted(), QuantLib.Time.Unadjusted(), QuantLib.Time.DateGenerationBackwards(), false,
                                        QuantLib.Time.USNYSECalendar())
fixing_days = 2
in_arrears = true
gearings = ones(length(fb_schedule.dates) - 1)
spreads = fill(0.001, length(fb_schedule.dates) - 1)
libor_3m = usd_libor_index(QuantLib.Time.TenorPeriod(Base.Dates.Month(3)), ts)
floating_bond = FloatingRateBond(settlement_days, face_amount, fb_schedule, libor_3m, fb_dc, conv, fixing_days, fb_issue_date, fb_bond_engine, in_arrears, 100.0, gearings, spreads, cap_vol=cap_vol)

# test
@test accrual_start_date(fixedrate_bond.cashflows.coupons[1]) == Date(2007, 5, 15)
@test accrual_end_date(fixedrate_bond.cashflows.coupons[1]) == Date(2007, 11, 15)
@test ref_period_start(fixedrate_bond.cashflows.coupons[1]) == Date(2007, 5, 15)
@test ref_period_end(fixedrate_bond.cashflows.coupons[1]) == Date(2007, 11, 15)
@test QuantLib.accrual_period(fixedrate_bond.cashflows.coupons[1]) == 0.5
@test QuantLib.get_latest_coupon(fixedrate_bond.cashflows) == fixedrate_bond.cashflows.coupons[end-1]

# test NPV calculations
@test npv(fixedrate_bond.cashflows, ts, ts.referenceDate, ts.referenceDate) == 94.67245111530346
@test npv(fixedrate_bond.cashflows, bbIR, true, ts.referenceDate, ts.referenceDate) == 94.71618240316593
@test QuantLib.npvbps(fixedrate_bond.cashflows, ts, ts.referenceDate, ts.referenceDate) == (94.67245111530346,0.07152580522996536)

# duration and yield
@test duration(ModifiedDuration(), fixedrate_bond.cashflows, bbIR, bbIR.dc, true, ts.referenceDate, ts.referenceDate) == 6.899980120600781
@test QuantLib.yield(fixedrate_bond.cashflows, 94.67245111530346, bbIR.dc, bbIR.comp, bbIR.freq, true, ts.referenceDate, ts.referenceDate, 1.0e-10, 100, 0.05) == 0.05506693320721388

# Other calculations
@test accrued_amount(fixedrate_bond.cashflows, ts.referenceDate) == 1.5407608695652275
@test previous_cashflow_date(fixedrate_bond.cashflows, ts.referenceDate) == Date(2008, 5, 15)
@test accrual_days(fixedrate_bond.cashflows, bbIR.dc, ts.referenceDate) == 126
@test next_cashflow(fixedrate_bond.cashflows, ts.referenceDate) == 3
@test next_coupon_rate(fixedrate_bond.cashflows, ts.referenceDate) == 0.045
@test has_occurred(fixedrate_bond.cashflows.coupons[end-1], ts.referenceDate) == false
@test maturity_date(fixedrate_bond.cashflows) == Date(2017, 5, 15)
@test amount(fixedrate_bond.cashflows.coupons[1]) == 2.2499999999999964
@test QuantLib.calc_rate(fixedrate_bond.cashflows.coupons[1]) == 0.045

# Floating coupon calculations
@test QuantLib.calc_rate(floating_bond.cashflows[12]) == 0.05443944335694899
@test amount(floating_bond.cashflows[12]) == 1.3912302191220296
@test QuantLib.index_fixing(floating_bond.cashflows[12]) == 0.05343944335694899
@test accrued_amount(floating_bond.cashflows[12], ts.referenceDate) == 0.8922019883499973
@test accrued_amount(floating_bond.cashflows[1], ts.referenceDate) == 0.0
