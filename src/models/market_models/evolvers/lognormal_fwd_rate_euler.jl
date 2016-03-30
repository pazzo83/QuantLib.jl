type LogNormalFwdRateEuler{M <: AbstractMarketModel, G <: BrownianGenerator} <: AbstractMarketModelEvolver
  marketModel::M
  numeraires::Vector{Int}
  initialStep::Int
  generator::G
  fixedDrifts::Vector{Vector{Float64}}
  numberOfRates::Int
  numberOfFactors::Int
  curveState::LMMCurveState
  currentStep::Int
  forwards::Vector{Float64}
  displacement::Vector{Float64}
  logForwards::Vector{Float64}
  initialLogForwards::Vector{Float64}
  drifts1::Vector{Float64}
  initialDrifts::Vector{Float64}
  brownians::Vector{Float64}
  correlatedBrownians::Vector{Float64}
  alive::Vector{Int}
  calculators::Vector{LMMDriftCalculator}
end

function LogNormalFwdRateEuler(marketModel::AbstractMarketModel, factory::BrownianGeneratorFactory, numeraires::Vector{Int}, initialStep::Int = 1)
  numberOfRates = marketModel.numberOfRates
  numberOfFactors = marketModel.numberOfFactors
  curveState = LMMCurveState(marketModel.evolution.rateTimes)
  forwards = copy(marketModel.initialRates)
  displacements = copy(marketModel.displacements)
  logForwards = Vector{Float64}(numberOfRates)
  initialLogForwards = Vector{Float64}(numberOfRates)
  drifts1 = Vector{Float64}(numberOfRates)
  initialDrifts = Vector{Float64}(numberOfRates)
  brownians = Vector{Float64}(numberOfFactors)
  correlatedBrownians = Vector{Float64}(numberOfRates)
  alive = marketModel.evolution.firstAliveRate

  check_compatibility(marketModel.evolution, numeraires)

  steps = number_of_steps(marketModel.evolution)

  generator = create(factory, numberOfFactors, steps - (initialStep-1))

  currentStep = initialStep

  calculators = Vector{LMMDriftCalculator}(steps)
  fixedDrifts = Vector{Vector{Float64}}(steps)

  for j = 1:steps
    A = marketModel.pseudoRoots[j]
    calculators[j] = LMMDriftCalculator(A, displacements, marketModel.evolution.rateTaus, numeraires[j], alive[j])
    fixed = Vector{Float64}(numberOfRates)

    for k in eachindex(fixed)
      # variance_ = dot(A[:, k], A[:, k])
      variance_ = vecdot(A[k, :], A[k, :])
      fixed[k] = -0.5 * variance_
    end

    fixedDrifts[j] = fixed
  end

  lognorm = LogNormalFwdRateEuler(marketModel, numeraires, initialStep, generator, fixedDrifts, numberOfRates, numberOfFactors, curveState, currentStep, forwards,
                              displacements, logForwards, initialLogForwards, drifts1, initialDrifts, brownians, correlatedBrownians, alive, calculators)

  set_forwards!(lognorm, marketModel.initialRates)
  return lognorm
end

function set_forwards!(lognorm::LogNormalFwdRateEuler, forwards::Vector{Float64})
  length(forwards) == lognorm.numberOfRates || error("mismatch between forwards and rate times")
  for i in eachindex(lognorm.initialLogForwards)
    lognorm.initialLogForwards[i] = log(forwards[i] + lognorm.displacement[i])
  end

  compute!(lognorm.calculators[lognorm.initialStep], forwards, lognorm.initialDrifts)

  return lognorm
end

function start_new_path!(lognorm::LogNormalFwdRateEuler)
  lognorm.currentStep = lognorm.initialStep
  lognorm.logForwards = copy(lognorm.initialLogForwards)

  return next_path!(lognorm.generator)
end

function advance_step!(lognorm::LogNormalFwdRateEuler)
  # we are going from T1 to T2

  # a) compute drifts D1 at T1
  if lognorm.currentStep > lognorm.initialStep
    compute!(lognorm.calculators[lognorm.currentStep], lognorm.forwards, lognorm.drifts1)
  else
    lognorm.drifts1 = copy(lognorm.initialDrifts)
  end

  # b) evolve forwards up to T2 using D1
  weight = next_step!(lognorm.generator, lognorm.brownians)
  A = copy(lognorm.marketModel.pseudoRoots[lognorm.currentStep])
  fixedDrift = copy(lognorm.fixedDrifts[lognorm.currentStep])

  alive = lognorm.alive[lognorm.currentStep]
  for i = alive:lognorm.numberOfRates
    lognorm.logForwards[i] += lognorm.drifts1[i] + fixedDrift[i]
    lognorm.logForwards[i] += vecdot(A[i, :], lognorm.brownians)
    lognorm.forwards[i] = exp(lognorm.logForwards[i]) - lognorm.displacement[i]
  end

  # same as PC evolver with two steps dropped

  # c) update curve state
  set_on_forward_rates!(lognorm.curveState, lognorm.forwards)
  lognorm.currentStep += 1

  return weight
end

current_state(lognorm::LogNormalFwdRateEuler) = lognorm.curveState
brownians_this_step(lognorm::LogNormalFwdRateEuler) = lognorm.brownians

# Clone #
clone(lognorm::LogNormalFwdRateEuler) = LogNormalFwdRateEuler(clone(lognorm.marketModel), copy(lognorm.numeraires), lognorm.initialStep, clone(lognorm.generator),
                                      deepcopy(lognorm.fixedDrifts), lognorm.numberOfRates, lognorm.numberOfFactors, clone(lognorm.curveState),
                                      lognorm.currentStep, copy(lognorm.forwards), copy(lognorm.displacement), copy(lognorm.logForwards), copy(lognorm.initialLogForwards),
                                      copy(lognorm.drifts1), copy(lognorm.initialDrifts), copy(lognorm.brownians), copy(lognorm.correlatedBrownians),
                                      copy(lognorm.alive), deepcopy(lognorm.calculators))
