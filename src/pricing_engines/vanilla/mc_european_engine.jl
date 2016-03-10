type MCEuropeanEngine{S <: AbstractBlackScholesProcess, I <: Integer} <: MCVanillaEngine
  process::S
  timeSteps::I
  timeStepsPerYear::I
  requiredSamples::I
  maxSamples::I
  requiredTolerance::Float64
  brownianBridge::Bool
  seed::I
  mcSimulation::MCSimulation
end

function MCEuropeanEngine{RSG <: AbstractRandomSequenceGenerator}(process::AbstractBlackScholesProcess; timeSteps::Int = -1, timeStepsPerYear::Int = -1, brownianBridge::Bool = false,
                          antitheticVariate::Bool = false, requiredSamples::Int = -1, requiredTolerance::Float64 = -1.0, maxSamples::Int = typemax(Int), seed::Int = 1, rsg::RSG = InverseCumulativeRSG(seed))
  # build mc sim
  mcSim = MCSimulation{RSG, SingleVariate}(antitheticVariate, false, rsg, SingleVariate())

  return MCEuropeanEngine(process, timeSteps, timeStepsPerYear, requiredSamples, maxSamples, requiredTolerance, brownianBridge, seed, mcSim)
end

type EuropeanPathPricer <: AbstractPathPricer
  payoff::PlainVanillaPayoff
  discount::Float64
end

EuropeanPathPricer(optionType::OptionType, strike::Float64, disc::Float64) = EuropeanPathPricer(PlainVanillaPayoff(optionType, strike), disc)

call(pricer::EuropeanPathPricer, path::Path) = pricer.payoff(path[end]) * pricer.discount

function time_grid(pe::MCVanillaEngine, opt::VanillaOption)
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

function path_generator(pe::MCVanillaEngine, opt::VanillaOption)
  dimensions = get_factors(pe.process)
  grid = time_grid(pe, opt)
  init_sequence_generator!(pe.mcSimulation.rsg, dimensions * (length(grid.times) - 1))

  return PathGenerator(pe.process, grid, pe.mcSimulation.rsg, pe.brownianBridge)
end

function path_pricer(pe::MCEuropeanEngine, opt::VanillaOption)
  payoff = opt.payoff
  process = pe.process

  return EuropeanPathPricer(payoff.optionType, payoff.strike, discount(process.riskFreeRate, time_grid(pe, opt)[end]))
end

function _calculate!(pe::MCVanillaEngine, opt::VanillaOption)
  _calculate!(pe.mcSimulation, pe, opt, pe.requiredTolerance, pe.requiredSamples, pe.maxSamples)
  opt.results.value = stats_mean(pe.mcSimulation.mcModel.sampleAccumulator)
  return pe, opt
end
