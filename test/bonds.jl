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
faceAmount = 100.0
fixed_schedule = QuantLib.Time.Schedule(Date(2007, 5, 15), Date(2017, 5, 15),
                       QuantLib.Time.TenorPeriod(QuantLib.Time.Semiannual()), QuantLib.Time.Unadjusted(),
                       QuantLib.Time.Unadjusted(), QuantLib.Time.DateGenerationBackwards(), false,
                       QuantLib.Time.USGovernmentBondCalendar())

pe = DiscountingBondEngine(ts)
fixedrate_bond = FixedRateBond(settlement_days, faceAmount, fixed_schedule, 0.045,
                        QuantLib.Time.ISMAActualActual(), QuantLib.Time.ModifiedFollowing(), 100.0,
                        Date(2007, 5, 15), fixed_schedule.cal, pe)


# generate a zero-coupon bond
zcb_pe = DiscountingBondEngine(ts)
zcb_cal = QuantLib.Time.USGovernmentBondCalendar()
zcb = ZeroCouponBond(3, zcb_cal, 100.0, Date(2013, 8, 15), QuantLib.Time.Following(), 116.92, Date(2003, 8, 15), zcb_pe)

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

# Zero coupon bond calculations
@test npv(zcb) == 89.54351524038537
@test clean_price(zcb) == 89.54351524038537
@test dirty_price(zcb) == 89.54351524038537
@test accrued_amount(zcb, settlement_date) == 0.0
@test QuantLib.yield(zcb, clean_price(zcb), QuantLib.Time.Actual360(), CompoundedCompounding(), QuantLib.Time.Annual(), settlement_date) == 0.05505323890596629

## Testing callable bonds ##
todaysDate = Date(2007, 10, 16)
set_eval_date!(settings, todaysDate)

ts = FlatForwardTermStructure(todaysDate, bbIR.rate, bbIR.dc, bbIR.comp, bbIR.freq)

# set up call schedule
callSchedule = CallabilitySchedule(24)
callPrice = 100.0
callDate = Date(2006, 9, 15)

for i in eachindex(callSchedule)
  nullCal = QuantLib.Time.NullCalendar()
  myPrice = Price(callPrice, CleanCall())

  callSchedule[i] = Callability(myPrice, Call(), callDate)
  callDate = QuantLib.Time.advance(Dates.Month(3), nullCal, callDate)
end

# set up callable bond
dated = Date(2004, 9, 16)
issue = dated
mat = Date(2012, 9, 15)

settlementDays = 3
bondCal = QuantLib.Time.USGovernmentBondCalendar()
coupon = 0.0465
freq = QuantLib.Time.Quarterly()
redemption = 100.0

bondDayCount = QuantLib.Time.ActualActualBond()
accrualConvention = QuantLib.Time.Unadjusted()
paymentConvention = QuantLib.Time.Unadjusted()

schedule = QuantLib.Time.Schedule(dated, mat, QuantLib.TenorPeriod(freq), accrualConvention, accrualConvention, QuantLib.Time.DateGenerationBackwards(), false, bondCal)

maxIter = 1000
accuracy = 1e-8
gridIntervals = 40
reversionParameter = 0.03

sigma = eps()
hw0 = HullWhite(ts, reversionParameter, sigma)
engine0 = TreeCallableFixedRateEngine(hw0, gridIntervals)

callableBond = CallableFixedRateBond(settlementDays, faceAmount, schedule, coupon, bondDayCount, paymentConvention, redemption, issue, callSchedule, engine0)

@test clean_price(callableBond) == 96.46794584230878
@test QuantLib.yield(callableBond, bondDayCount, CompoundedCompounding(), freq, accuracy, maxIter) == 0.054753332138061515
