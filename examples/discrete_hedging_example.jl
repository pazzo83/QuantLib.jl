using QuantLib

## TYPES ##
type ReplicationError{OT <: OptionType}
  maturity::Float64
  payoff::PlainVanillaPayoff{OT}
  s0::Float64
  sigma::Float64
  r::Float64
  vega::Float64
end

function ReplicationError{OT <: OptionType}(optionType::OT, maturity::Float64, strike::Float64, s0::Float64, sigma::Float64, r::Float64)
  rDiscount = exp(-r * maturity)
  qDiscount = 1.0
  forward = s0 * rDiscount / rDiscount
  stdDev = sqrt(sigma * sigma * maturity)
  payoff = PlainVanillaPayoff(optionType, strike)

  black = BlackCalculator(payoff, forward, stdDev, rDiscount)

  vega_ = vega(black, maturity)

  return ReplicationError{OT}(maturity, payoff, s0, sigma, r, vega_)
end

type ReplicationPathPricer{OT <: OptionType} <: QuantLib.AbstractPathPricer
  optionType::OT
  strike::Float64
  r::Float64
  maturity::Float64
  sigma::Float64
end

## Type Methods ##

function (rpp::ReplicationPathPricer)(path::Path)
  n = length(path)
  # discrete hedging interval
  dt = rpp.maturity / n
  # we assume the stock pays no dividend
  stockDividendYield = 0.0
  # starting...
  t = 0.0
  stock = path[1] # stock at t=0
  money_account = 0.0 # money_account at t=0

  ### Initial Deal ###
  # Option fair price (black-scholes) at t=0
  rDiscount = exp(-rpp.r * rpp.maturity)
  qDiscount = exp(-stockDividendYield * rpp.maturity)
  forward = stock * qDiscount / rDiscount
  stdDev = sqrt(rpp.sigma * rpp.sigma * rpp.maturity)

  payoff = PlainVanillaPayoff(rpp.optionType, rpp.strike)
  black = BlackCalculator(payoff, forward, stdDev, rDiscount)
  # sell the option, cash in its premium
  money_account += value(black)
  # compute delta
  delta_ = delta(black, stock)
  # delta-hedge the option buying stock
  stockAmount = delta_
  money_account -= stockAmount * stock

  ## Heding during option life ##
  for step = 1:n-1
    # time flows
    t += dt
    # accruing on the money account
    money_account *= exp(rpp.r * dt)

    # stock growth
    stock = path[step + 1]

    # recalculate the option value at the current stock value, and the current time to maturity
    rDiscount = exp(-rpp.r * (rpp.maturity - t))
    qDiscount = exp(-stockDividendYield * (rpp.maturity - t))
    forward = stock * qDiscount / rDiscount
    stdDev = sqrt(rpp.sigma * rpp.sigma * (rpp.maturity - t))
    black = BlackCalculator(payoff, forward, stdDev, rDiscount)
    # recalculate delta
    delta_ = delta(black, stock)
    # re-hedging
    money_account -= (delta_ - stockAmount) * stock
    stockAmount = delta_
  end

  ### Option expiration ###
  # last accrual on my money account
  money_account *= exp(rpp.r * dt)
  # last stock growth
  stock = path[n]
  # the hedger delivers the option payoff to the option holder
  optionPayoff = PlainVanillaPayoff(rpp.optionType, rpp.strike)(stock)
  money_account -= optionPayoff
  # and unwinds the hedge selling his stock position
  money_account += stockAmount * stock
  # final PnL
  return money_account
end

function compute(rp::ReplicationError, nTimeSteps::Int, nSamples::Int)
  cal = QuantLib.Time.TargetCalendar()
  tdy = Dates.today()

  dc = QuantLib.Time.Actual365()

  stateVar = Quote(rp.s0)
  riskFreeRate = FlatForwardTermStructure(tdy, rp.r, dc)
  dividendYield = FlatForwardTermStructure(tdy, 0.0, dc)
  vol = BlackConstantVol(tdy, cal, rp.sigma, dc)
  diffuse = BlackScholesMertonProcess(stateVar, riskFreeRate, dividendYield, vol)

  # Black scholes equation rules the path generator
  # At each step of the log of the stock will ahve drift and sigma^2 variance
  rsg = QuantLib.Math.InverseCumulativeRSG(0, nTimeSteps)
  brownianBridge = false
  myPathGenerator = PathGenerator(diffuse, rp.maturity, nTimeSteps, rsg, brownianBridge)

  myPathPricer = ReplicationPathPricer{typeof(rp.payoff.optionType)}(rp.payoff.optionType, rp.payoff.strike, rp.r, rp.maturity, rp.sigma)
  statsAccumulator = QuantLib.Math.gen_RiskStatistics()

  mcSim = MonteCarloModel(myPathGenerator, myPathPricer, statsAccumulator, false)

  add_samples!(mcSim, nSamples)

  PLMean = QuantLib.Math.stats_mean(mcSim.sampleAccumulator)
  PLStdDev = QuantLib.Math.stats_std_deviation(mcSim.sampleAccumulator)
  PLSkew = QuantLib.Math.stats_skewness(mcSim.sampleAccumulator)
  PLKurtosis = QuantLib.Math.stats_kurtosis(mcSim.sampleAccumulator)

  println("P&L Mean: ", PLMean)
  println("P&L Standard Deviation: ", PLStdDev)
  println("P&L Skew: ", PLSkew)
  println("P&L Kurtosis: ", PLKurtosis)
end

function main()
  maturity = 1.0 / 12.0 # 1 month
  strike = 100.0
  underlying = 100.0
  vol = 0.20 # 20%
  riskFreeRate = 0.05 # 5%
  rp = ReplicationError(Call(), maturity, strike, underlying, vol, riskFreeRate)

  scenarios = 50000
  hedgesNum = 21
  compute(rp, hedgesNum, scenarios)
  hedgesNum = 84
  compute(rp, hedgesNum, scenarios)
end
