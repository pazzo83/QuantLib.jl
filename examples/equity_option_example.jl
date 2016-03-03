# Equity Option example
# various pricing methods

using QuantLib

function main()
  cal = QuantLib.Time.TargetCalendar()
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
  dc = QuantLib.Time.Actual365()

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

  europeanOption = update_pricing_engine(europeanOption, hestonPE)
  println(npv(europeanOption))

  # bates process
  batesProcess = BatesProcess(flatTermStructure, flatDividendTS, underlyingH, vol * vol, 1.0, vol * vol,
                              0.001, 0.0, 1e-14, 1e-14, 1e-14)

  batesModel = BatesModel(batesProcess)

  batesPE = BatesEngine(batesModel)

  europeanOption = update_pricing_engine(europeanOption, batesPE)
  println(npv(europeanOption))

  # American Option
  baroneAdesiPE = BaroneAdesiWhaleyApproximationEngine(bsmProcess)
  americanOption = VanillaOption(payoff, americanExercise, baroneAdesiPE)
  println(npv(americanOption))

  # Bjerksund Stensland engine
  bjerksundStenslandPE = BjerksundStenslandApproximationEngine(bsmProcess)
  americanOption = update_pricing_engine(americanOption, bjerksundStenslandPE)
  println(npv(americanOption))

  # Integral Engine
  integralPE = IntegralEngine(bsmProcess)
  europeanOption = update_pricing_engine(europeanOption, integralPE)
  println(npv(europeanOption))

  # FD Engines
  timeSteps = 801
  fdEuroPE = FDEuropeanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
  europeanOption = update_pricing_engine(europeanOption, fdEuroPE)
  println(npv(europeanOption))

  # Bermudan Option
  fdBermudanPE = FDBermudanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
  bermudanOption = VanillaOption(payoff, bermudanExercise, fdBermudanPE)
  println(npv(bermudanOption))

  fdAmericanPE = FDAmericanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
  americanOption = update_pricing_engine(americanOption, fdAmericanPE)
end
