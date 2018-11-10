using QuantLib

function theVegaBumps(factorwiseBumping::Bool, marketModel::AbstractMarketModel, doCaps::Bool)
  multiplierCutoff = 50.0
  projectionTolerance = 1e-4
  numberRates = marketModel.numberOfRates
  caps = VolatilityBumpInstrumentJacobianCap[]

  if doCaps
    capStrike = marketModel.initialRates[1]
    for i = 1:numberRates - 1
      nextCap = VolatilityBumpInstrumentJacobianCap(i, i, capStrike)
      push!(caps, nextCap)
    end
  end

  swaptions = Vector{VolatilityBumpInstrumentJacobianSwaption}(undef, numberRates)

  for i in eachindex(swaptions)
    swaptions[i] = VolatilityBumpInstrumentJacobianSwaption(i, numberRates)
  end

  possibleBumps = VegaBumpCollection(marketModel, factorwiseBumping)

  bumpFinder = OrthogonalizedBumpFinder(possibleBumps, swaptions, caps, multiplierCutoff, projectionTolerance)

  theBumps = Vector{Vector{Matrix{Float64}}}()
  get_vega_bumps!(bumpFinder, theBumps)
end

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
  basisCoefficients = Vector{Vector{Float64}}()

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
  clonedEvolver = QuantLib.clone(evolver)

  t1 = time()

  collect_node_data!(evolver, inverseFloater, basisSystem, nullRebate, control, trainingPaths, collectedData)

  t2 = time()

  # calcualte the exercise stategy's coefficients
  generic_longstaff_schwartz_regression!(collectedData, basisCoefficients)

  # turn the coefficients into an exercise strategy
  exerciseStrategy = LongstaffSchwartzExerciseStrategy(basisSystem, basisCoefficients, evolution, numeraires, nullRebate, control)

  # callable receiver swap
  callableProduct = CallSpecifiedMultiProduct(inverseFloater, exerciseStrategy, ExerciseAdapter(nullRebate))

  allProducts = MarketModelComposite()
  add_product!(allProducts, inverseFloater)
  add_product!(allProducts, callableProduct)
  finalize!(allProducts)

  accounter = AccountingEngine(clonedEvolver, QuantLib.clone(allProducts), initialNumeraireValue)

  stats = QuantLib.Math.GenericSequenceStats()
  multiple_path_values!(accounter, stats, paths)

  t3 = time()

  means = QuantLib.Math.stats_mean(stats)

  println("Means: ", means)
  println("Time to build strategy: ", t2 - t1, " seconds")
  println("Time to price: ", t3 - t2, " seconds")

  # vegas

  # do it twice once with factorwise bumping, one without
  pathsToVegas = vegaPaths

  for i = 0:3
    allowFactorwiseBumping = i % 2 > 0
    doCaps = round(Int, i / 2) > 0

    evolverEuler = LogNormalFwdRateEuler(marketModel, generatorFactory, numeraires)

    pathwiseInverseFloater = MarketModelPathwiseInverseFloater(rateTimes, accruals, accruals, fixedStrikes, fixedMultipliers, floatingSpreads, paymentTimes, payer)

    pathwiseInverseFloaterClone = QuantLib.clone(pathwiseInverseFloater)

    callableProductPathwise = CallSpecifiedPathwiseMultiProduct(pathwiseInverseFloaterClone, exerciseStrategy)

    theBumps = theVegaBumps(allowFactorwiseBumping, marketModel, doCaps)

    accountingEngineVegas = PathwiseVegasOuterAccountingEngine(QuantLib.clone(evolverEuler), QuantLib.clone(callableProductPathwise), marketModel, theBumps, initialNumeraireValue)

    values = Vector{Float64}()
    errors = Vector{Float64}()

    multiple_path_values!(accountingEngineVegas, values, errors, pathsToVegas)

    println("vega output")
    println("factorwise bumping: ", allowFactorwiseBumping)
    println("doCaps: ", doCaps)

    r = 1
    println("price estimate: ", values[r])
    r += 1

    for i = 1:numberOfRates
      println("Delta ", i, ", ", values[r], ", ", errors[r])
      r += 1
    end

    totalVega = 0.0

    for t = r:length(values)
      println("vega, ", t - 1 - numberOfRates, ", ", values[t], ", ", errors[t])
      totalVega += values[t]
    end

    println("total vega: ", totalVega)
  end

  ## TODO Upper Bound
end

function main()
  for i = 5:9
    InverseFloater(i / 100.0)
  end
end
