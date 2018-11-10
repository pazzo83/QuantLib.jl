# Equity Option example
# various pricing methods

using QuantLib
using Dates

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
  println("BlackScholes (Euro):                   ", npv(europeanOption))

  # heston process
  hestonProcess = HestonProcess(flatTermStructure, flatDividendTS, underlyingH, vol * vol, 1.0, vol * vol, 0.001, 0.0)
  hestonModel = HestonModel(hestonProcess)
  hestonPE = AnalyticHestonEngine(hestonModel)

  europeanOption = update_pricing_engine(europeanOption, hestonPE)
  println("Heston semi-analytic (Euro):           ", npv(europeanOption))

  # bates process
  batesProcess = BatesProcess(flatTermStructure, flatDividendTS, underlyingH, vol * vol, 1.0, vol * vol,
                              0.001, 0.0, 1e-14, 1e-14, 1e-14)
  batesModel = BatesModel(batesProcess)
  batesPE = BatesEngine(batesModel)

  europeanOption = update_pricing_engine(europeanOption, batesPE)
  println("Bates semi-analytic (Euro):            ", npv(europeanOption))

  # American Option
  baroneAdesiPE = BaroneAdesiWhaleyApproximationEngine(bsmProcess)
  americanOption = VanillaOption(payoff, americanExercise, baroneAdesiPE)
  println("Barone-Adesi/Whaley (Amer):            ", npv(americanOption))

  # Bjerksund Stensland engine
  bjerksundStenslandPE = BjerksundStenslandApproximationEngine(bsmProcess)
  americanOption = update_pricing_engine(americanOption, bjerksundStenslandPE)
  println("Bjerksund Stensland (Amer):            ", npv(americanOption))

  # Integral Engine
  integralPE = IntegralEngine(bsmProcess)
  europeanOption = update_pricing_engine(europeanOption, integralPE)
  println("Integral (Euro):                       ", npv(europeanOption))

  # FD Engines
  timeSteps = 801
  fdEuroPE = FDEuropeanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
  europeanOption = update_pricing_engine(europeanOption, fdEuroPE)
  println("FD (Euro):                             ", npv(europeanOption))

  # Bermudan Option
  fdBermudanPE = FDBermudanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
  bermudanOption = VanillaOption(payoff, bermudanExercise, fdBermudanPE)
  println("FD (Berm):                             ", npv(bermudanOption))

  fdAmericanPE = FDAmericanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
  americanOption = update_pricing_engine(americanOption, fdAmericanPE)
  println("FD (Amer):                             ", npv(americanOption))

  # Binomial Tree engines
  # Jarrow Rudd
  btJarrowRuddEuro = BinomialVanillaEngine(bsmProcess, timeSteps, JarrowRudd)
  europeanOption = update_pricing_engine(europeanOption, btJarrowRuddEuro)
  println("Binomial Jarrow-Rudd (Euro):           ", npv(europeanOption))

  btJarrowRuddBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, JarrowRudd)
  bermudanOption = update_pricing_engine(bermudanOption, btJarrowRuddBermudan)
  println("Binomial Jarrow-Rudd (Berm):           ", npv(bermudanOption))

  btJarrowRuddAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, JarrowRudd)
  americanOption = update_pricing_engine(americanOption, btJarrowRuddAmerican)
  println("Binomial Jarrow-Rudd (Amer):           ", npv(americanOption))

  # CoxRossRubinstein
  btCoxRubinsteinEuro = BinomialVanillaEngine(bsmProcess, timeSteps, CoxRossRubinstein)
  europeanOption = update_pricing_engine(europeanOption, btCoxRubinsteinEuro)
  println("Binomial Cox-Ross-Rubinstein (Euro):   ", npv(europeanOption))

  btCoxRubinsteinBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, CoxRossRubinstein)
  bermudanOption = update_pricing_engine(bermudanOption, btCoxRubinsteinBermudan)
  println("Binomial Cox-Ross-Rubinstein (Berm):   ", npv(bermudanOption))

  btCoxRubinsteinAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, CoxRossRubinstein)
  americanOption = update_pricing_engine(americanOption, btCoxRubinsteinAmerican)
  println("Binomial Cox-Ross-Rubinstein (Amer):   ", npv(americanOption))

  # AdditiveEQP
  btAdditiveEQPEuro = BinomialVanillaEngine(bsmProcess, timeSteps, AdditiveEQP)
  europeanOption = update_pricing_engine(europeanOption, btAdditiveEQPEuro)
  println("Additive Equiprobs (Euro):             ", npv(europeanOption))

  btAdditiveEQPBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, AdditiveEQP)
  bermudanOption = update_pricing_engine(bermudanOption, btAdditiveEQPBermudan)
  println("Additive Equiprobs (Berm):             ", npv(bermudanOption))

  btAdditiveEQPAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, AdditiveEQP)
  americanOption = update_pricing_engine(americanOption, btAdditiveEQPAmerican)
  println("Additive Equiprobs (Amer):             ", npv(americanOption))

  # Trigeorgis
  btTrigeorgisEuro = BinomialVanillaEngine(bsmProcess, timeSteps, Trigeorgis)
  europeanOption = update_pricing_engine(europeanOption, btTrigeorgisEuro)
  println("Binomial Trigeorgis (Euro):            ", npv(europeanOption))

  btTrigeorgisBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, Trigeorgis)
  bermudanOption = update_pricing_engine(bermudanOption, btTrigeorgisBermudan)
  println("Binomial Trigeorgis (Berm):            ", npv(bermudanOption))

  btTrigeorgisAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, Trigeorgis)
  americanOption = update_pricing_engine(americanOption, btTrigeorgisAmerican)
  println("Binomial Trigeorgis (Amer):            ", npv(americanOption))

  # Tian
  btTianEuro = BinomialVanillaEngine(bsmProcess, timeSteps, Tian)
  europeanOption = update_pricing_engine(europeanOption, btTianEuro)
  println("Binomial Tian (Euro):                  ", npv(europeanOption))

  btTianBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, Tian)
  bermudanOption = update_pricing_engine(bermudanOption, btTianBermudan)
  println("Binomial Tian (Berm):                  ", npv(bermudanOption))

  btTianAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, Tian)
  americanOption = update_pricing_engine(americanOption, btTianAmerican)
  println("Binomial Tian (Amer):                  ", npv(americanOption))

  # LeisenReimer
  btLeisenReimerEuro = BinomialVanillaEngine(bsmProcess, timeSteps, LeisenReimer)
  europeanOption = update_pricing_engine(europeanOption, btLeisenReimerEuro)
  println("Binomial Leisen-Reimer (Euro):         ", npv(europeanOption))

  btLeisenReimerBermudan = BinomialVanillaEngine(bsmProcess, timeSteps, LeisenReimer)
  bermudanOption = update_pricing_engine(bermudanOption, btLeisenReimerBermudan)
  println("Binomial Leisen-Reimer (Berm):         ", npv(bermudanOption))

  btLeisenReimerAmerican = BinomialVanillaEngine(bsmProcess, timeSteps, LeisenReimer)
  americanOption = update_pricing_engine(americanOption, btLeisenReimerAmerican)
  println("Binomial Leisen-Reimer (Amer):         ", npv(americanOption))

  # Joshi4
  btJoshi4Euro = BinomialVanillaEngine(bsmProcess, timeSteps, Joshi4)
  europeanOption = update_pricing_engine(europeanOption, btJoshi4Euro)
  println("Binomial Joshi (Euro):                 ", npv(europeanOption))

  btJoshi4Bermudan = BinomialVanillaEngine(bsmProcess, timeSteps, Joshi4)
  bermudanOption = update_pricing_engine(bermudanOption, btJoshi4Bermudan)
  println("Binomial Joshi (Berm):                 ", npv(bermudanOption))

  btJoshi4American = BinomialVanillaEngine(bsmProcess, timeSteps, Joshi4)
  americanOption = update_pricing_engine(americanOption, btJoshi4American)
  println("Binomial Joshi (Amer):                 ", npv(americanOption))

  ## MONTE CARLO ##
  timeSteps = 1
  mcSeed = 42
  mcengine1 = MCEuropeanEngine(bsmProcess; timeSteps = timeSteps, requiredTolerance = 0.02, seed = mcSeed)
  europeanOption = update_pricing_engine(europeanOption, mcengine1)
  println("MC-crude (Euro):                       ", npv(europeanOption))

  nSamples = 32768
  mcengine2 = MCEuropeanEngine(bsmProcess; timeSteps = timeSteps, requiredSamples = nSamples)
  europeanOption = update_pricing_engine(europeanOption, mcengine2)
  println("MC-Sobol (Euro):                       ", npv(europeanOption))

  mcengine3 = MCAmericanEngine(bsmProcess; timeSteps = 100, antitheticVariate=true, requiredTolerance=0.02, seed=mcSeed, nCalibrationSamples=4096)
  americanOption = update_pricing_engine(americanOption, mcengine3)
  println("MC-Longstaff Schwartz (Amer):          ", npv(americanOption))
  # americanOption
end
