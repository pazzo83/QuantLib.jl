# Equity Option example
# various pricing methods

using JQuantLib

function main()
  cal = JQuantLib.Time.TargetCalendar()
  todaysDate = Date(1998, 5, 15)
  settlementDate = Date(1998, 5, 17)

  set_eval_date!(settings, todaysDate)

  optionType = Put()
  underlying = 36.0
  strike = 40.0
  dividendYield = 0.00
  riskFreeRate = 0.06
  vol = 0.20
  mat = Date(1999, 5, 17)
  dc = JQuantLib.Time.Actual365()

  exerciseDates = Date[settlementDate + Dates.Month(3 * i) for i = 1:4]

  europeanExercise = EuropeanExercise(mat)
  bermudanExercise = BermudanExercise(exerciseDates)
  americanExercise = AmericanExercise(settlementDate, mat)

  underlyingH = Quote(underlying)

  flatTermStructure = FlatForwardTermStructure(settlementDate, cal, Quote(riskFreeRate), dc)
  flatDividendTS = FlatForwardTermStructure(settlementDate, cal, Quote(dividendYield), dc)
  flatVolTS = BlackConstantVol(settlementDate, cal, vol, dc)

  payoff = PlainVanillaPayoff(optionType, strike)

  bsmProcess = BlackScholesMertonProcess(underlyingH, flatTermStructure, flatDividendTS, flatVolTS)

  # options
  bsPE = AnalyticEuropeanEngine(bsmProcess) # black scholes PE
  europeanOption = VanillaOption(payoff, europeanExercise, bsPE)

  # black scholes for European
  println(npv(europeanOption))

  # heston process
  hestonProcess = HestonProcess(flatTermStructure, flatDividendTS, underlyingH, vol * vol, 1.0, vol * vol, 0.001, 0.0)

  hestonModel = HestonModel(hestonProcess)

  hestonPE = AnalyticHestonEngine(hestonModel)
end
