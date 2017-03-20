type MCEuropeanEngine{S <: AbstractBlackScholesProcess, RSG <: AbstractRandomSequenceGenerator} <: MCVanillaEngine{S, RSG}
  process::S
  timeSteps::Int
  timeStepsPerYear::Int
  requiredSamples::Int
  maxSamples::Int
  requiredTolerance::Float64
  brownianBridge::Bool
  seed::Int
  # mcSimulation::MCSimulation{RSG, T}
  antitheticVariate::Bool
  rsg::RSG
end

function MCEuropeanEngine{S <: AbstractBlackScholesProcess, RSG <: AbstractRandomSequenceGenerator}(process::S; timeSteps::Int = -1, timeStepsPerYear::Int = -1, brownianBridge::Bool = false,
                          antitheticVariate::Bool = false, requiredSamples::Int = -1, requiredTolerance::Float64 = -1.0, maxSamples::Int = typemax(Int), seed::Int = 1, rsg::RSG = InverseCumulativeRSG(seed))
  # build mc sim
  # mcSim = MCSimulation{RSG, SingleVariate}(antitheticVariate, false, rsg, SingleVariate())

  return MCEuropeanEngine{S, RSG}(process, timeSteps, timeStepsPerYear, requiredSamples, maxSamples, requiredTolerance, brownianBridge, seed, antitheticVariate, rsg)
end

type EuropeanPathPricer{OT <: OptionType} <: AbstractPathPricer
  payoff::PlainVanillaPayoff{OT}
  discount::Float64
end

EuropeanPathPricer{OT <: OptionType}(optionType::OT, strike::Float64, disc::Float64) = EuropeanPathPricer{OT}(PlainVanillaPayoff(optionType, strike), disc)

(pricer::EuropeanPathPricer)(path::Path) = pricer.payoff(path[end]) * pricer.discount

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
  init_sequence_generator!(pe.rsg, dimensions * (length(grid.times) - 1))

  return PathGenerator(pe.process, grid, pe.rsg, pe.brownianBridge)
end

function path_pricer(pe::MCEuropeanEngine, opt::VanillaOption)
  payoff = opt.payoff
  process = pe.process

  return EuropeanPathPricer(payoff.optionType, payoff.strike, discount(process.riskFreeRate, time_grid(pe, opt)[end]))
end

function _calculate!(pe::MCVanillaEngine, opt::VanillaOption)
  mcSim = MCSimulation(pe, false, opt, SingleVariate())
  _calculate!(mcSim, pe.requiredTolerance, pe.requiredSamples, pe.maxSamples)
  opt.results.value = stats_mean(mcSim.mcModel.sampleAccumulator)
  return pe, opt
end
