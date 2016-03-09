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
  mcSim
end

function MCAmericanEngine{RSG <: AbstractRandomSequenceGenerator}(process::AbstractBlackScholesProcess; timeSteps::Int = -1, timeStepsPerYear::Int = -1, brownianBridge::Bool = false,
                          antitheticVariate::Bool = false, requiredSamples::Int = -1, requiredTolerance::Float64 = -1.0, maxSamples::Int = typemax(Int), seed::Int = 0, rsg::RSG = InverseCumulativeRSG(seed),
                          nCalibrationSamples::Int = 2048, polynomOrder::Int = 2, polynomType::LsmBasisSystemPolynomType = Monomial())
  # build mc sim
  mcSim = MCSimulation{RSG, SingleVariate}(antitheticVariate, false, rsg, gen_RiskStatistics(), SingleVariate())

  return MCAmericanEngine(process, timeSteps, timeStepsPerYear, requiredSamples, maxSamples, requiredTolerance, brownianBridge, seed, nCalibrationSamples, polynomOrder, polynomType, mcSim)
end
