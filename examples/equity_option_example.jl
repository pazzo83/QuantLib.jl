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
  println(npv(americanOption))

  # Binomial Tree engines
  # Jarrow Rudd
  btJarrowRuddEuro = BinomialVanillaEngine(bsmProcess, timeSteps, JarrowRudd)
  europeanOption = update_pricing_engine(europeanOption, btJarrowRuddEuro)
  println(npv(europeanOption))

  btJarrowRuddBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, JarrowRudd)
  bermudanOption = update_pricing_engine(bermudanOption, btJarrowRuddBermudan)
  println(npv(bermudanOption))

  btJarrowRuddAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, JarrowRudd)
  americanOption = update_pricing_engine(americanOption, btJarrowRuddAmerican)
  println(npv(americanOption))

  # CoxRossRubinstein
  btCoxRubinsteinEuro = BinomialVanillaEngine(bsmProcess, timeSteps, CoxRossRubinstein)
  europeanOption = update_pricing_engine(europeanOption, btCoxRubinsteinEuro)
  println(npv(europeanOption))

  btCoxRubinsteinBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, CoxRossRubinstein)
  bermudanOption = update_pricing_engine(bermudanOption, btCoxRubinsteinBermudan)
  println(npv(bermudanOption))

  btCoxRubinsteinAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, CoxRossRubinstein)
  americanOption = update_pricing_engine(americanOption, btCoxRubinsteinAmerican)
  println(npv(americanOption))

  # AdditiveEQP
  btAdditiveEQPEuro = BinomialVanillaEngine(bsmProcess, timeSteps, AdditiveEQP)
  europeanOption = update_pricing_engine(europeanOption, btAdditiveEQPEuro)
  println(npv(europeanOption))

  btAdditiveEQPBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, AdditiveEQP)
  bermudanOption = update_pricing_engine(bermudanOption, btAdditiveEQPBermudan)
  println(npv(bermudanOption))

  btAdditiveEQPAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, AdditiveEQP)
  americanOption = update_pricing_engine(americanOption, btAdditiveEQPAmerican)
  println(npv(americanOption))

  # Trigeorgis
  btTrigeorgisEuro = BinomialVanillaEngine(bsmProcess, timeSteps, Trigeorgis)
  europeanOption = update_pricing_engine(europeanOption, btTrigeorgisEuro)
  println(npv(europeanOption))

  btTrigeorgisBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, Trigeorgis)
  bermudanOption = update_pricing_engine(bermudanOption, btTrigeorgisBermudan)
  println(npv(bermudanOption))

  btTrigeorgisAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, Trigeorgis)
  americanOption = update_pricing_engine(americanOption, btTrigeorgisAmerican)
  println(npv(americanOption))

  # Tian
  btTianEuro = BinomialVanillaEngine(bsmProcess, timeSteps, Tian)
  europeanOption = update_pricing_engine(europeanOption, btTianEuro)
  println(npv(europeanOption))

  btTianBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, Tian)
  bermudanOption = update_pricing_engine(bermudanOption, btTianBermudan)
  println(npv(bermudanOption))

  btTianAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, Tian)
  americanOption = update_pricing_engine(americanOption, btTianAmerican)
  println(npv(americanOption))

  # LeisenReimer
  btLeisenReimerEuro = BinomialVanillaEngine(bsmProcess, timeSteps, LeisenReimer)
  europeanOption = update_pricing_engine(europeanOption, btLeisenReimerEuro)
  println(npv(europeanOption))

  btLeisenReimerBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, LeisenReimer)
  bermudanOption = update_pricing_engine(bermudanOption, btLeisenReimerBermudan)
  println(npv(bermudanOption))

  btLeisenReimerAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, LeisenReimer)
  americanOption = update_pricing_engine(americanOption, btLeisenReimerAmerican)
  println(npv(americanOption))

  # Joshi4
  btJoshi4Euro = BinomialVanillaEngine(bsmProcess, timeSteps, Joshi4)
  europeanOption = update_pricing_engine(europeanOption, btJoshi4Euro)
  println(npv(europeanOption))

  btJoshi4Bermudan = BinomialVanillaEngine(bsmProcess, timeSteps, Joshi4)
  bermudanOption = update_pricing_engine(bermudanOption, btJoshi4Bermudan)
  println(npv(bermudanOption))

  btJoshi4American = BinomialVanillaEngine(bsmProcess, timeSteps, Joshi4)
  americanOption = update_pricing_engine(americanOption, btJoshi4American)
  println(npv(americanOption))

  ## MONTE CARLO ##
  timeSteps = 1
  mcSeed = 42
  mcengine1 = MCEuropeanEngine(bsmProcess; timeSteps = timeSteps, requiredTolerance = 0.02, seed = mcSeed)
  europeanOption = update_pricing_engine(europeanOption, mcengine1)
  println(npv(europeanOption))

  nSamples = 32768
  mcengine2 = MCEuropeanEngine(bsmProcess; timeSteps = timeSteps, requiredSamples = nSamples)
  europeanOption = update_pricing_engine(europeanOption, mcengine2)
  println(npv(europeanOption))

  # mcengine3 = MCAmericanEngine(bsmProcess; timeSteps = 100, antitheticVariate=true, requiredTolerance=0.02, seed=mcSeed, nCalibrationSamples=4096)
  # americanOption = update_pricing_engine(americanOption, mcengine3)
  # println(npv(americanOption))
  # europeanOption
end
