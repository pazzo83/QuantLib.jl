using QuantLib

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
    swapT = underlying.swapT == QuantLib.Receiver() ? "Receiver" : "Payer"
    println("$expiry    $endDate      $nominal       $rate   $swapT      $vol")
  end
end

function print_model_calibration(basket::Vector{SwaptionHelper}, volatilities::Vector{Float64})
  println("Expiry   Model Sigma  ModelPx  MktPx     ModelVol  MktVol")
  println("=========================================================")
  for i in eachindex(basket)
    helper = basket[i]
    expiry = helper.swaption.exercise.dates[1]
    modelSig = volatilities[i]
    modelPx = model_value!(helper)
    mktPx = helper.calibCommon.marketValue
    modelVol = implied_volatility!(helper, modelPx, 1e-6, 1000, 0.0, 2.0)
    marketVol = helper.volatility.value

    println(@sprintf("%s  %.6f  %.6f  %.6f  %.2f     %.2f", expiry, modelSig, modelPx, mktPx, modelVol, marketVol))
  end
end

function main()
  refDate = Dates.Date(2014, 4, 30)
  set_eval_date!(settings, refDate)

  forward6mLevel = 0.025
  oisLevel = 0.02

  forward6mQuote = Quote(forward6mLevel)
  oisQuote = Quote(oisLevel)

  yts6m = FlatForwardTermStructure(0, QuantLib.Time.TargetCalendar(), forward6mQuote, QuantLib.Time.Actual365())
  ytsOis = FlatForwardTermStructure(0, QuantLib.Time.TargetCalendar(), oisQuote, QuantLib.Time.Actual365())

  euribor6m = euribor_index(QuantLib.Time.TenorPeriod(Dates.Month(6)), yts6m)

  volLevel = 0.2
  volQuote = Quote(volLevel)

  swaptionVol = ConstantSwaptionVolatility(0, QuantLib.Time.TargetCalendar(), QuantLib.Time.ModifiedFollowing(), volQuote, QuantLib.Time.Actual365())

  strike = 0.04

  effectiveDate = QuantLib.Time.advance(Dates.Day(2), QuantLib.Time.TargetCalendar(), refDate)
  maturityDate = QuantLib.Time.advance(Dates.Year(10), QuantLib.Time.TargetCalendar(), effectiveDate)

  schedBizConvention = QuantLib.Time.ModifiedFollowing()

  fixedSchedule = QuantLib.Time.Schedule(effectiveDate, maturityDate, QuantLib.Time.TenorPeriod(Dates.Year(1)), schedBizConvention, schedBizConvention,
                  QuantLib.Time.DateGenerationForwards(), false, QuantLib.Time.TargetCalendar())
  floatingSchedule = QuantLib.Time.Schedule(effectiveDate, maturityDate, QuantLib.Time.TenorPeriod(Dates.Month(6)), schedBizConvention, schedBizConvention,
                  QuantLib.Time.DateGenerationForwards(), false, QuantLib.Time.TargetCalendar())

  exerciseDates = Vector{Date}(9)

  for i = 2:10
    exerciseDates[i - 1] = QuantLib.Time.advance(Dates.Day(-2), QuantLib.Time.TargetCalendar(), fixedSchedule.dates[i])
  end

  stepDates = exerciseDates[1:end - 1]
  sigmas = fill(0.01, length(exerciseDates))
  reversion = 0.01

  gsr = GSR(yts6m, stepDates, sigmas, reversion)

  swaptionEngine = Gaussian1DSwaptionEngine(gsr, 64, 7.0, true, false, ytsOis)
  nonstandardSwaptionEngine = Gaussian1DNonstandardSwaptionEngine(gsr, 64, 7.0, true, false, Quote(-1.0), ytsOis)

  underlying = NonstandardSwap(VanillaSwap(Payer(), 1.0, fixedSchedule, strike, QuantLib.Time.BondThirty360(), euribor6m, 0.0, floatingSchedule, QuantLib.Time.Actual360()))

  exercise = BermudanExercise(exerciseDates)

  swaption = NonstandardSwaption(underlying, exercise, nonstandardSwaptionEngine)

  swapBase = EuriborSwapIsdaFixA(QuantLib.Time.TenorPeriod(Dates.Year(10)), yts6m, ytsOis)

  basket = calibration_basket(swaption, nonstandardSwaptionEngine, swapBase, swaptionVol, NaiveBasketType())

  print_calibration_basket(basket)

  for i in eachindex(basket)
    basket[i] = update_pricing_engine(basket[i], swaptionEngine)
  end

  method = QuantLib.Math.LevenbergMarquardt()

  ec = QuantLib.Math.EndCriteria(1000, 10, 1.0e-8, 1.0e-8, 1.0e-8)

  calibrate_volatilities_iterative!(gsr, basket, method, ec)

  print_model_calibration(basket, get_volatilities(gsr))

  println(npv(swaption))
end
