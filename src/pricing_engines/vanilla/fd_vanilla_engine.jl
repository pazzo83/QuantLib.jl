type FDEuropeanEngine{B <: AbstractBlackScholesProcess, I <: Integer} <: AbstractFDVanillaEngine
  process::B
  timeSteps::I
  gridPoints::I
  timeDependent::Bool
  exerciseDate::Date
  finiteDifferenceOperator::TridiagonalOperator
  intrinsicValues::SampledCurve
  BCs::Vector{BoundaryCondition}
  sMin::Float64
  center::Float64
  sMax::Float64
  prices::SampledCurve
  fdEvolverFunc::Function
end

function FDEuropeanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100,
                          gridPoints::Int = 100, timeDependent::Bool = false)
  prices = SampledCurve(gridPoints)
  finiteDifferenceOperator, intrinsicValues, BCs = gen_fd_vanilla_engine_params(gridPoints)
  sMin = center = sMax = 0.0
  exerciseDate = Date()

  return FDEuropeanEngine(process, timeSteps, gridPoints, timeDependent, exerciseDate, finiteDifferenceOperator,
                          intrinsicValues, BCs, sMin, center, sMax, prices, fdEvolverFunc)
end

function gen_fd_vanilla_engine_params(gridPoints::Int)
  finiteDifferenceOperator = TridiagonalOperator()
  intrinsicValues = SampledCurve(gridPoints)
  BCs = Vector{BoundaryCondition}(2)

  return finiteDifferenceOperator, intrinsicValues, BCs
end

get_residual_time(pe::AbstractFDVanillaEngine) = get_time(pe.process, pe.exerciseDate)

function safe_grid_points(gridPoints::Int, residualTime::Float64)
  minGridPoints = 10
  minGridPointsPerYear = 2
  return max(gridPoints, residualTime > 1.0 ? round(Int, floor(minGridPoints + (residualTime - 1.0) * minGridPointsPerYear)) : minGridPoints)
end

function set_grid_limits!(pe::AbstractFDVanillaEngine, center::Float64, t::Float64)
  center > 0.0 || error("negative or null underlying given")
  t > 0.0 || error("negative or zero residual time")

  pe.center = center
  newGridPoints = safe_grid_points(pe.gridPoints, t)
  if newGridPoints > JQuantLib.Math.get_size(pe.intrinsicValues)
    pe.intrinsicValues = SampledCurve(newGridPoints)
  end

  volSqrtTime = sqrt(black_variance(pe.process.blackVolatility, t, pe.center))

  # the prefactor fine-tunes performance at small volatilities
  prefactor = 1.0 + 0.02 / volSqrtTime
  minMaxFactor = exp(4.0 * prefactor * volSqrtTime)

  pe.sMin = pe.center / minMaxFactor
  pe.sMax = pe.center * minMaxFactor

  return pe
end

function ensure_strike_in_grid!(pe::AbstractFDVanillaEngine, opt::VanillaOption)
  striked_payoff = opt.payoff
  requiredGridValue = striked_payoff.strike
  safetyZoneFactor = 1.1

  if pe.sMin > requiredGridValue / safetyZoneFactor
    pe.sMin = requiredGridValue / safetyZoneFactor
    # enforce central placement of underlying
    pe.sMax = pe.center / (pe.sMin / pe.center)
  end

  if pe.sMax < requiredGridValue * safetyZoneFactor
    pe.sMax = requiredGridValue * safetyZoneFactor
    # enforce central placement of underlying
    pe.sMin = pe.center / (pe.sMax / pe.center)
  end

  return pe
end

function set_grid_limits!(pe::AbstractFDVanillaEngine, opt::VanillaOption)
  set_grid_limits!(pe, state_variable(pe.process).value, get_residual_time(pe))
  ensure_strike_in_grid!(pe, opt)

  return pe
end

function initialize_initial_condition!(pe::AbstractFDVanillaEngine, opt::VanillaOption)
  set_log_grid!(pe.intrinsicValues, pe.sMin, pe.sMax)
  sample!(pe.intrinsicValues, opt.payoff)

  return pe
end

function initialize_boundary_conditions!(pe::AbstractFDVanillaEngine)
  pe.BCs[1] = NeumannBC(pe.intrinsicValues.values[2] - pe.intrinsicValues.values[1], LowerSide())
  pe.BCs[2] = NeumannBC(pe.intrinsicValues.values[end] - pe.intrinsicValues.values[end - 1], UpperSide())

  return pe
end

initialize_operator!(pe::AbstractFDVanillaEngine) = pe.finiteDifferenceOperator = get_operator(pe.process, pe.intrinsicValues.grid, get_residual_time(pe), pe.timeDependent)

function _calculate!(pe::FDEuropeanEngine, opt::EuropeanOption)
  pe.exerciseDate = opt.exercise.dates[end]
  payoff = opt.payoff

  set_grid_limits!(pe, opt)
  initialize_initial_condition!(pe, opt)
  initialize_operator!(pe)
  initialize_boundary_conditions!(pe)

  model = FiniteDifferenceModel(pe.finiteDifferenceOperator, pe.BCs, pe.fdEvolverFunc)
  return pe, opt
end
