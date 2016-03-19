type LogNormalFwdRatePc{M <: AbstractMarketModel, G <: BrownianGenerator} <: AbstractMarketModelEvolver
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
  drifts2::Vector{Float64}
  initialDrifts::Vector{Float64}
  brownians::Vector{Float64}
  correlatedBrownians::Vector{Float64}
  alive::Vector{Int}
  calculators::Vector{LMMDriftCalculator}
end

function LogNormalFwdRatePc(marketModel::AbstractMarketModel, factory::BrownianGeneratorFactory, numeraires::Vector{Int}, initialStep::Int = 1)
  numberOfRates = marketModel.numberOfRates
  numberOfFactors = marketModel.numberOfFactors
  curveState = LMMCurveState(marketModel.evolution.rateTimes)
  forwards = marketModel.initialRates
  displacements = marketModel.displacements
  logForwards = Vector{Float64}(numberOfRates)
  initialLogForwards = Vector{Float64}(numberOfRates)
  drifts1 = Vector{Float64}(numberOfRates)
  drifts2 = Vector{Float64}(numberOfRates)
  initialDrifts = Vector{Float64}(numberOfRates)
  brownians = Vector{Float64}(numberOfFactors)
  correlatedBrownians = Vector{Float64}(numberOfRates)
  alive = marketModel.evolution.firstAliveRate

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

  lognorm = LogNormalFwdRatePc(marketModel, numeraires, initialStep, generator, fixedDrifts, numberOfRates, numberOfFactors, curveState, currentStep, forwards,
                              displacements, logForwards, initialLogForwards, drifts1, drifts2, initialDrifts, brownians, correlatedBrownians, alive, calculators)

  set_forwards!(lognorm, marketModel.initialRates)
  return lognorm
end

current_state(lognorm::LogNormalFwdRatePc) = lognorm.curveState

function set_forwards!(lognorm::LogNormalFwdRatePc, forwards::Vector{Float64})
  length(forwards) == lognorm.numberOfRates || error("mismatch between forwards and rate times")
  for i in eachindex(lognorm.initialLogForwards)
    lognorm.initialLogForwards[i] = log(forwards[i] + lognorm.displacement[i])
  end

  compute!(lognorm.calculators[lognorm.initialStep], forwards, lognorm.initialDrifts)

  return lognorm
end

function start_new_path!(lognorm::LogNormalFwdRatePc)
  lognorm.currentStep = lognorm.initialStep
  lognorm.logForwards = copy(lognorm.initialLogForwards)

  return next_path!(lognorm.generator)
end

function advance_step!(lognorm::LogNormalFwdRatePc)
  # we're going from T1 to T2

  # a) compute drifts D1 at T1
  if lognorm.currentStep > lognorm.initialStep
    compute!(lognorm.calculators[lognorm.currentStep], lognorm.forwards, lognorm.drifts1)
  else
    lognorm.drifts1 = copy(lognorm.initialDrifts)
  end

  # b) evolve forwards up to T2 using D1
  weight = next_step!(lognorm.generator, lognorm.brownians)
  A = lognorm.marketModel.pseudoRoots[lognorm.currentStep]
  fixedDrift = lognorm.fixedDrifts[lognorm.currentStep]

  alive = lognorm.alive[lognorm.currentStep]
  for i = alive:lognorm.numberOfRates
    lognorm.logForwards[i] += lognorm.drifts1[i] + fixedDrift[i]
    lognorm.logForwards[i] += vecdot(A[i, :], lognorm.brownians)
    lognorm.forwards[i] = exp(lognorm.logForwards[i]) - lognorm.displacement[i]
  end

  # c) recompute drifts D2 using the predicted forwards
  compute!(lognorm.calculators[lognorm.currentStep], lognorm.forwards, lognorm.drifts2)

  # d) correct forwards using both drifts
  for i = alive:lognorm.numberOfRates
    lognorm.logForwards[i] += (lognorm.drifts2[i] - lognorm.drifts1[i]) / 2.0
    lognorm.forwards[i] = exp(lognorm.logForwards[i]) - lognorm.displacement[i]
  end

  # e) update curve state
  set_on_forward_rates!(lognorm.curveState, lognorm.forwards)
  lognorm.currentStep += 1

  return weight
end
