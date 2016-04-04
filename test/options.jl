## This also tests some vanilla pricing engines ##
using Base.Test
using QuantLib

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

# test a European Option
bsPE = AnalyticEuropeanEngine(bsmProcess) # black scholes PE
europeanOption = VanillaOption(payoff, europeanExercise, bsPE)

# American Option
baroneAdesiPE = BaroneAdesiWhaleyApproximationEngine(bsmProcess)
americanOption = VanillaOption(payoff, americanExercise, baroneAdesiPE)

# Bermudan Option
timeSteps = 801
fdBermudanPE = FDBermudanEngine(bsmProcess, CrankNelson, timeSteps, timeSteps - 1)
bermudanOption = VanillaOption(payoff, bermudanExercise, fdBermudanPE)

@test npv(europeanOption) == 3.8443077915968398
@test npv(americanOption) == 4.459627613776478
@test npv(bermudanOption) == 4.3608068625088645
