type MCAmericanEngine{S <: AbstractBlackScholesProcess, I <: Integer, P <: LsmBasisSystemPolynomType} <: MCLongstaffSchwartzEngine
  process::S
  timeSteps::I
  timeStepsPerYear::I
  requiredSamples::I
  maxSamples::I
  requiredTolerance::Float64
  brownianBridge::Bool
  seed::I
  nCalibrationSamples::I
  polynomOrder::I
  polynomType::P
  mcSimulation::MCSimulation
end

function MCAmericanEngine{RSG <: AbstractRandomSequenceGenerator}(process::AbstractBlackScholesProcess; timeSteps::Int = -1, timeStepsPerYear::Int = -1, brownianBridge::Bool = false,
                          antitheticVariate::Bool = false, requiredSamples::Int = -1, requiredTolerance::Float64 = -1.0, maxSamples::Int = typemax(Int), seed::Int = 0, rsg::RSG = InverseCumulativeRSG(seed),
                          nCalibrationSamples::Int = 2048, polynomOrder::Int = 2, polynomType::LsmBasisSystemPolynomType = Monomial())
  # build mc sim
  mcSim = MCSimulation{RSG, SingleVariate}(antitheticVariate, false, rsg, SingleVariate())

  return MCAmericanEngine(process, timeSteps, timeStepsPerYear, requiredSamples, maxSamples, requiredTolerance, brownianBridge, seed, nCalibrationSamples, polynomOrder, polynomType, mcSim)
end

type AmericanPathPricer{P <: StrikedTypePayoff, T, T1} <: EarlyExercisePathPricer
  scalingValue::Float64
  payoff::P
  v::Vector{Union{T, T1}}
end

function AmericanPathPricer{P <: StrikedTypePayoff, L <: LsmBasisSystemPolynomType}(payoff::P, polynomOrder::Int, polynomType::L)
  T = get_type(polynomType)
  v = Vector{Union{T, P}}(polynomOrder + 2)
  path_basis_system!(polynomType, polynomOrder, v)
  v[end] = payoff

  scalingVal = 1.0 / payoff.strike

  return AmericanPathPricer{P, T, P}(scalingVal, payoff, v)
end

basis_system(p::AmericanPathPricer) = p.v
get_state(p::AmericanPathPricer, path::Path, t::Int) = path[t] * p.scalingValue
get_payoff(p::AmericanPathPricer, st::Float64) = p.payoff(st / p.scalingValue)
call(p::AmericanPathPricer, path::Path, t::Int) = get_payoff(p, get_state(path, t))

type LongstaffSchwartzPathPricer{E <: EarlyExercisePathPricer, I <: Integer, T} <: AbstractPathPricer
  calibrationPhase::Bool
  pathPricer::E
  coeff::Vector{Vector{Float64}}
  dF::Vector{Float64}
  paths::Vector{Path}
  len::I
  exerciseProbability::NonWeightedStatistics
  v::Vector{T}
end

function LongstaffSchwartzPathPricer(tg::TimeGrid, ep::EarlyExercisePathPricer, yts::YieldTermStructure)
  v = basis_system(ep)
  len = length(tg.times)
  coeff = Vector{Vector{Float64}}(len - 1)
  dF = zeros(len - 1)
  paths = Vector{Path}()

  for i in eachindex(dF)
    dF[i] = discount(yts, tg.times[i+1]) / discount(yts, tg.times[i])
  end

  return LongstaffSchwartzPathPricer(true, ep, coeff, dF, paths, len, NonWeightedStatistics(), v)
end

function call(lpp::LongstaffSchwartzPathPricer, path::Path)
  if lpp.calibrationPhase
    push!(lpp.paths, path)
    return 0.0
  end

  price = lpp.pathPricer(path, len)

  # initialize with exercise on last date
  exercised = price != 0.0

  for i = len-1:-1:2
    price *= lpp.dF[i]
    exercise = lpp.pathPricer(path, i)
    if exercise > 0.0
      regValue = get_state(lpp.pathPricer, path, i)

      continuationValue = 0.0
      for l in eachindex(lpp.v)
        continuationValue += lpp.coeff[i-1][l] * lpp.v[l](regValue)
      end

      if continuationValue < exercise
        price = exercise
        exercised = true
      end
    end
  end
  add_sample!(lpp.exerciseProbability, exercised ? 1.0 : 0.0)

  return price * lpp.dF[1]
end

function time_grid(pe::MCLongstaffSchwartzEngine, opt::VanillaOption)
  lastExerciseDate = opt.exercise.dates[end]
  t = get_time(pe.process, lastExerciseDate)
  if pe.timeSteps != -1
    return TimeGrid(t, pe.timeSteps)
  elseif pe.timeStepsPerYear != -1
    steps = round(Int, floor(pe.timeStepsPerYear * t))
    return TimeGrid(t, max(steps, 1))
  else
    error("time steps not specified")
  end
end

function lsm_path_pricer(pe::MCAmericanEngine, opt::AmericanOption)
  process = pe.process
  exercise = opt.exercise

  early_exercise_path_pricer = AmericanPathPricer(opt.payoff, pe.polynomOrder, pe.polynomType)

  return LongstaffSchwartzPathPricer(time_grid(pe, opt), early_exercise_path_pricer, process.riskFreeRate)
end

function path_generator(pe::MCLongstaffSchwartzEngine, opt::VanillaOption)
  dimensions = get_factors(pe.process)
  grid = time_grid(pe, opt)
  init_sequence_generator!(pe.mcSimulation.rsg, dimensions * (length(grid.times) - 1))

  return PathGenerator(pe.process, grid, pe.mcSimulation.rsg, pe.brownianBridge)
end

function _calculate!(pe::MCLongstaffSchwartzEngine, opt::VanillaOption)
  path_pricer = lsm_path_pricer(pe, opt)
  pe.mcSimulation.mcModel = MonteCarloModel(path_generator(pe, opt), path_pricer, gen_RiskStatistics(), pe.mcSimulation.antitheticVariate)
  add_samples!(pe.mcSimulation.mcModel, pe.nCalibrationSamples, 1)

  return pe, opt
end
