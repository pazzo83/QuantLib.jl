using QuantLib

function main()
  todaysDate = Date(2007, 10, 16)
  set_eval_date!(settings, todaysDate)

  bbCurveRate = 0.055
  bbDayCount = QuantLib.Time.ActualActualBond()
  bbIR = InterestRate(bbCurveRate, bbDayCount, CompoundedCompounding(), QuantLib.Time.Semiannual())
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
  coupon = [0.0465]
  freq = QuantLib.Time.Quarterly()
  redemption = 100.0
  faceAmount = 100.0

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

  println(clean_price(callableBond))
  println(QuantLib.yield(callableBond, bondDayCount, CompoundedCompounding(), freq, accuracy, maxIter))
end
