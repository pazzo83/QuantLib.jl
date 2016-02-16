using JQuantLib

function print_calibration_basket(basket::Vector{SwaptionHelper})
  println("Expiry         Maturity       Nominal        Rate               Pay/Rec    Market ivol")
  println("======================================================================================")
  for i in eachindex(basket)
    helper = basket[i]
    underlying = underlying_swap!(helper)
    endDate = underlying.fixedSchedule.dates[end]
    nominal = underlying.nominal
    vol = helper.volatility.value
    rate = underlying.fixedRate
    expiry = helper.swaption.exercise.dates[1]
    swapT = underlying.swapT == JQuantLib.Receiver() ? "Receiver" : "Payer"
    println("$expiry    $endDate      $nominal       $rate   $swapT      $vol")
  end
end

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

  exerciseDates = Vector{Date}(9)

  for i = 2:10
    exerciseDates[i - 1] = JQuantLib.Time.advance(Dates.Day(-2), JQuantLib.Time.TargetCalendar(), fixedSchedule.dates[i])
  end

  stepDates = exerciseDates[1:end - 1]
  sigmas = fill(0.01, length(exerciseDates))
  reversion = 0.01

  gsr = GSR(yts6m, stepDates, sigmas, reversion)

  swaptionEngine = Gaussian1DSwaptionEngine(gsr, 64, 7.0, true, false, ytsOis)
  nonstandardSwaptionEngine = Gaussian1DNonstandardSwaptionEngine(gsr, 64, 7.0, true, false, Quote(-1.0), ytsOis)

  underlying = NonstandardSwap(VanillaSwap(Payer(), 1.0, fixedSchedule, strike, JQuantLib.Time.BondThirty360(), euribor6m, 0.0, floatingSchedule, JQuantLib.Time.Actual360()))

  exercise = BermudanExercise(exerciseDates)

  swaption = NonstandardSwaption(underlying, exercise, nonstandardSwaptionEngine)

  swapBase = EuriborSwapIsdaFixA(JQuantLib.Time.TenorPeriod(Dates.Year(10)), yts6m, ytsOis)

  basket = calibration_basket(swaption, nonstandardSwaptionEngine, swapBase, swaptionVol, NaiveBasketType())

  print_calibration_basket(basket)

  for i in eachindex(basket)
    basket[i] = update_pricing_engine(basket[i], swaptionEngine)
  end

  method = JQuantLib.Math.LevenbergMarquardt()

  ec = JQuantLib.Math.EndCriteria(1000, 10, 1.0e-8, 1.0e-8, 1.0e-8)

  calibrate_volatilities_iterative!(gsr, basket, method, ec)
end
