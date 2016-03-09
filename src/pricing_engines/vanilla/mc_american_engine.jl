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

type LongstaffSchwartzPathPricer{E <: EarlyExercisePathPricer, T} <: AbstractPathPricer
  calibrationPhase::Bool
  pathPricer::E
  coeff::Vector{Float64}
  dF::Vector{Float64}
  paths::Vector{Path}
  v::Vector{T}
end

function LongstaffSchwartzPathPricer(tg::TimeGrid, ep::EarlyExercisePathPricer, yts::YieldTermStructure)
  v = basis_system(ep)
  coeff = zeros(length(tg.times) - 1)
  dF = zeros(length(tg.times) - 1)
  paths = Vector{Path}()

  return LongstaffSchwartzPathPricer(true, ep, coeff, dF, paths, v)
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

function _calculate!(pe::MCLongstaffSchwartzEngine, opt::VanillaOption)
  path_pricer = lsm_path_pricer(pe, opt)

  return pe, opt
end
