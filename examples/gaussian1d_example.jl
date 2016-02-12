using JQuantLib

function main()
  refDate = Dates.Date(2014, 4, 30)
  set_eval_date!(settings, refDate)

  forward6mLevel = 0.025
  oisLevel = 0.02

  forward6mQuote = Quote(forward6mLevel)
  oisQuote = Quote(oisLevel)

  yts6m = FlatForwardTermStructure(0, JQuantLib.Time.TargetCalendar(), forward6mQuote, JQuantLib.Time.Actual365())
  ytsOis = FlatForwardTermStructure(0, JQuantLib.Time.TargetCalendar(), oisQuote, JQuantLib.Time.Actual365())

  euribor6m = euribor_index(JQuantLib.Time.TenorPeriod(Dates.Month(6)), yts6m)

  volLevel = 0.2
  volQuote = Quote(volLevel)

  swaptionVol = ConstantSwaptionVolatility(0, JQuantLib.Time.TargetCalendar(), JQuantLib.Time.ModifiedFollowing(), volQuote, JQuantLib.Time.Actual365())

  strike = 0.04

  effectiveDate = JQuantLib.Time.advance(Dates.Day(2), JQuantLib.Time.TargetCalendar(), refDate)
  maturityDate = JQuantLib.Time.advance(Dates.Year(10), JQuantLib.Time.TargetCalendar(), effectiveDate)

  schedBizConvention = JQuantLib.Time.ModifiedFollowing()

  fixedSchedule = JQuantLib.Time.Schedule(effectiveDate, maturityDate, JQuantLib.Time.TenorPeriod(Dates.Year(1)), schedBizConvention, schedBizConvention,
                  JQuantLib.Time.DateGenerationForwards(), false, JQuantLib.Time.TargetCalendar())
  floatingSchedule = JQuantLib.Time.Schedule(effectiveDate, maturityDate, JQuantLib.Time.TenorPeriod(Dates.Month(6)), schedBizConvention, schedBizConvention,
                  JQuantLib.Time.DateGenerationForwards(), false, JQuantLib.Time.TargetCalendar())

  underlying = NonstandardSwap(VanillaSwap(Payer(), 1.0, fixedSchedule, strike, JQuantLib.Time.BondThirty360(), euribor6m, 0.0, floatingSchedule, JQuantLib.Time.Actual360()))

  exerciseDates = Vector{Date}(9)

  for i = 2:10
    exerciseDates[i - 1] = JQuantLib.Time.advance(Dates.Day(-2), JQuantLib.Time.TargetCalendar(), fixedSchedule.dates[i])
  end

  exercise = BermudanExercise(exerciseDates)

  swaption = NonstandardSwaption(underlying, exercise)

  return swaption
end
