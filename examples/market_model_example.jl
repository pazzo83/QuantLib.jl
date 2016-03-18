using QuantLib

function InverseFloater(rateLevel::Float64)
  numberOfRates = 20
  accrual = 0.5
  firstTime = 0.5

  strike = 0.15
  fixedMultiplier = 2.0
  floatingSpread = 0.0
  payer = true

  rateTimes = Float64[firstTime + (i-1) * accrual for i = 1:numberOfRates + 1]

  accruals = fill(accrual, numberOfRates)
  fixedStrikes = fill(strike, numberOfRates)
  floatingSpreads = fill(floatingSpread, numberOfRates)
  fixedMultipliers = fill(fixedMultiplier, numberOfRates)

  paymentTimes = Float64[firstTime + i * accrual for i = 1:numberOfRates]

  inverseFloater = MultiStepInverseFloater(rateTimes, accruals, accruals, fixedStrikes, fixedMultipliers, floatingSpreads, paymentTimes, payer)

  # exercise schedule, we can exercise on any rate time except the last one
  exerciseTimes = rateTimes[1:end-1]

  # naive exercise strategy, exercise above a trigger level
  trigger = 0.05
  swapTriggers = fill(trigger, length(exerciseTimes))
  naifStrategy = SwapRateTrigger(rateTimes, swapTriggers, exerciseTimes)

  # Longstaff-Schwartz exercise strategy
  collectedData = Vector{Vector{NodeData}}()
  basisCoefficients = Vector{Float64}()

  # control that does nothing, need it because some control is required
  control = NothingExerciseValue(rateTimes)

  basisSystem = SwapForwardBasisSystem(rateTimes, exerciseTimes)

  # rebate that does nothing, need it because some rebate is expected
  # when you break a swap, nothing happens
  nullRebate = NothingExerciseValue(rateTimes)

  dummyProduct = CallSpecifiedMultiProduct(inverseFloater, naifStrategy, ExerciseAdapter(nullRebate))

  evolution = dummyProduct.evolution

  # params for models
  seed = 12332
  trainingPaths = 65536
  paths = 65536
  vegaPaths = 16384

  println("inverse floater")
  println("fixed strikes: ", strike)
  println("number rates: ", numberOfRates)
  println("training paths: ", trainingPaths)
  println("paths: ", paths)
  println("vega paths: ", vegaPaths)

  # set up a callibration, this would typically be done by using a callibrator
  println("rate level: ", rateLevel)

  initialNumeraireValue = 0.95
  volLevel = 0.11
  beta = 0.2
  gamma = 1.0
  numberOfFactors = min(5, numberOfRates)
  displacementLevel = 0.02

  # set up vectors
  initialRates = fill(rateLevel, numberOfRates)
  volatilities = fill(volLevel, numberOfRates)
  displacements = fill(displacementLevel, numberOfRates)

  correlations = ExponentialForwardCorrelation(rateTimes, volLevel, beta, gamma)

  calibration = FlatVol(volatilities, QuantLib.clone(correlations), evolution, numberOfFactors, initialRates, displacements)
  marketModel = QuantLib.clone(calibration)

  # use a factory
  generatorFactory = SobolBrownianGeneratorFactory(SobolDiagonalOrdering(), seed)

  numeraires = money_market_measure(evolution)
  evolver = LogNormalFwdRatePc(marketModel, generatorFactory, numeraires)

  collect_node_data!(evolver, inverseFloater, basisSystem, nullRebate, control, trainingPaths, collectedData)
end

function main()
  InverseFloater(5 / 100.0)
end
