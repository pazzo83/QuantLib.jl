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

function MCEuropeanEngine(process::AbstractBlackScholesProcess, timeSteps::Int, timeStepsPerYear::Int, brownianBridge::Bool, antitheticVariate::Bool,
                          requiredSamples::Int, requiredTolerance::Float64, maxSamples::Int, seed::Int, rng::AbstractRandomNumberGenerator = PseudoRandom())
  # build mc sim
  mcSim = MCSimulation(antitheticVariate, false, rng, Statistics(), SingleVariate())

  return MCEuropeanEngine(process, timeSteps, timeStepsPerYear, requiredSamples, maxSamples, requiredTolerance, brownianBridge, seed, mcSim)
end
