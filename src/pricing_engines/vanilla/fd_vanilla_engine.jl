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

type FDBermudanEngine{B <: AbstractBlackScholesProcess, I <: Integer} <: FDMultiPeriodEngine
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
  stoppingTimes::Vector{Float64}
  timeStepPerPeriod::I
  fdEvolverFunc::Function
end

function FDBermudanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100,
                          gridPoints::Int = 100, timeDependent::Bool = false)
  prices = SampledCurve(gridPoints)
  finiteDifferenceOperator, intrinsicValues, BCs = gen_fd_vanilla_engine_params(gridPoints)
  sMin = center = sMax = 0.0
  exerciseDate = Date()
  stoppingTimes = Vector{Float64}()
  timeStepPerPeriod = timeSteps

  return FDBermudanEngine(process, timeSteps, gridPoints, timeDependent, exerciseDate, finiteDifferenceOperator,
                          intrinsicValues, BCs, sMin, center, sMax, prices, stoppingTimes, timeStepPerPeriod, fdEvolverFunc)
end

type FDAmericanEngine{B <: AbstractBlackScholesProcess, I <: Integer} <: FDStepConditionEngine
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
  controlPrices::SampledCurve
  controlBCs::Vector{BoundaryCondition}
  controlOperator::TridiagonalOperator
  fdEvolverFunc::Function
end

function FDAmericanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100,
                          gridPoints::Int = 100, timeDependent::Bool = false)
  # build engine
  prices = SampledCurve()
  controlPrices = SampledCurve(gridPoints)
  finiteDifferenceOperator, intrinsicValues, BCs = gen_fd_vanilla_engine_params(gridPoints)
  sMin = center = sMax = 0.0
  exerciseDate = Date()
  controlBCs = Vector{BoundaryCondition}(2)
  controlOperator = TridiagonalOperator()

  return FDAmericanEngine(process, timeSteps, gridPoints, timeDependent, exerciseDate, finiteDifferenceOperator,
                          intrinsicValues, BCs, sMin, center, sMax, prices, controlPrices, controlBCs,
                          controlOperator, fdEvolverFunc)
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
  if newGridPoints > QuantLib.Math.get_size(pe.intrinsicValues)
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
  pe.BCs[1] = build_NeumannBC(pe.intrinsicValues.values[2] - pe.intrinsicValues.values[1], LowerSide())
  pe.BCs[2] = build_NeumannBC(pe.intrinsicValues.values[end] - pe.intrinsicValues.values[end - 1], UpperSide())

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

  pe.prices = copy(pe.intrinsicValues)

  rollback!(model, pe.prices.values, get_residual_time(pe), 0.0, pe.timeSteps)
  opt.results.value = value_at_center(pe.prices)
  opt.results.delta = first_derivative_at_center(pe.prices)
  opt.results.gamma = second_derivative_at_center(pe.prices)
  opt.results.theta = black_scholes_theta(pe.process, opt.results.value, opt.results.delta, opt.results.gamma)
  return pe, opt
end

function setup_args!(pe::FDMultiPeriodEngine, opt::VanillaOption)
  n = length(opt.exercise.dates)

  pe.stoppingTimes = zeros(n)
  for i in eachindex(pe.stoppingTimes)
    pe.stoppingTimes[i] = get_time(pe.process, opt.exercise.dates[i])
  end

  return pe
end

get_dividend_time(pe::FDMultiPeriodEngine, idx::Int) = pe.stoppingTimes[idx]

function execute_intermediate_step!(pe::FDBermudanEngine, ::Int)
  for i in eachindex(pe.intrinsicValues.grid)
    pe.prices.values[i] = max(pe.prices.values[i], pe.intrinsicValues.values[i])
  end

  return pe
end

function _calculate!(pe::FDMultiPeriodEngine, opt::VanillaOption)
  setup_args!(pe, opt)
  pe.exerciseDate = opt.exercise.dates[end]
  payoff = opt.payoff

  dateNumber = length(pe.stoppingTimes)
  lastDateIsResTime = false
  firstIndex = 0
  lastIndex = dateNumber
  firstDateIsZero = false
  firstNonZeroDate = get_residual_time(pe)

  dateTolerance = 1e-6

  resid_time = get_residual_time(pe)

  if dateNumber > 0
    get_dividend_time(pe, 1) >= 0 || error("first date cannot be negative")
    if get_dividend_time(pe, 1) < resid_time * dateTolerance
      firstDateIsZero = true
      firstIndex = 1
      if dateNumber >= 2
        firstNonZeroDate = get_dividend_time(pe, 2)
      end
    end

    if abs(get_dividend_time(pe, lastIndex) - resid_time) < dateTolerance
      lastDateIsResTime = true
      lastIndex = dateNumber - 1
    end

    if ~firstDateIsZero
      firstNonZeroDate = get_dividend_time(pe, 1)
    end

    if dateNumber >= 2
      issorted(pe.stoppingTimes) || error("dates must be in increasing order")
    end
  end

  dt = resid_time / (pe.timeStepPerPeriod * (dateNumber + 1.0))

  # ensure that dt is always smaller than the first non-zero date
  if firstNonZeroDate <= dt
    dt = firstNonZeroDate / 2.0
  end

  set_grid_limits!(pe, opt)
  initialize_initial_condition!(pe, opt)
  initialize_operator!(pe)
  initialize_boundary_conditions!(pe)
  model = FiniteDifferenceModel(pe.finiteDifferenceOperator, pe.BCs, pe.fdEvolverFunc)

  pe.prices = copy(pe.intrinsicValues)

  if lastDateIsResTime
    execute_intermediate_step!(pe, dateNumber - 1)
  end

  j = lastIndex
  while true
    if j == dateNumber
      beginDate = get_residual_time(pe)
    else
      beginDate = get_dividend_time(pe, j+1)
    end

    if j >= 1
      endDate = get_dividend_time(pe, j)
    else
      endDate = dt
    end

    rollback!(model, pe.prices.values, beginDate, endDate, pe.timeStepPerPeriod)

    if j >= 1
      execute_intermediate_step!(pe, j)
    end

    j -= 1
    if j < firstIndex
      break
    end
  end

  rollback!(model, pe.prices.values, dt, 0.0, 1)

  if firstDateIsZero
    execute_intermediate_step!(pe, 0)
  end

  opt.results.value = value_at_center(pe.prices)
  opt.results.delta = first_derivative_at_center(pe.prices)
  opt.results.gamma = second_derivative_at_center(pe.prices)

  return pe, opt
end

initialize_step_condition(pe::FDAmericanEngine) = build_AmericanStepCondition(pe.intrinsicValues.values)

function _calculate!(pe::FDStepConditionEngine, opt::VanillaOption)
  pe.exerciseDate = opt.exercise.dates[end]
  payoff = opt.payoff

  set_grid_limits!(pe, opt)
  initialize_initial_condition!(pe, opt)
  initialize_operator!(pe)
  initialize_boundary_conditions!(pe)
  stepCond = initialize_step_condition(pe)

  pe.prices = copy(pe.intrinsicValues)
  pe.controlPrices = copy(pe.intrinsicValues)
  pe.controlOperator = copy(pe.finiteDifferenceOperator)
  pe.controlBCs[1] = pe.BCs[1]
  pe.controlBCs[2] = pe.BCs[2]

  operatorSet = TridiagonalOperator[pe.finiteDifferenceOperator, pe.controlOperator]
  arraySet = Vector{Float64}[x for x in (pe.prices.values, pe.controlPrices.values)]
  bcSet = hcat(pe.BCs, pe.controlBCs)
  conditionSet = StepConditionSet(StepCondition[stepCond, FdmNullStepCondition()])

  model = FiniteDifferenceModel(operatorSet, bcSet, pe.fdEvolverFunc)

  rollback!(model, arraySet, get_residual_time(pe), 0.0, pe.timeSteps, conditionSet)

  pe.prices.values = arraySet[1]
  pe.controlPrices.values = arraySet[2]

  variance = black_variance(pe.process.blackVolatility, opt.exercise.dates[end], payoff.strike)

  dividendDiscount = discount(pe.process.dividendYield, opt.exercise.dates[end])

  riskFreeDiscount = discount(pe.process.riskFreeRate, opt.exercise.dates[end])

  spot = state_variable(pe.process).value
  forwardPrice = spot *  dividendDiscount / riskFreeDiscount

  black = BlackCalculator(payoff, forwardPrice, sqrt(variance), riskFreeDiscount)

  opt.results.value = value_at_center(pe.prices) - value_at_center(pe.controlPrices) + value(black)
  opt.results.delta = first_derivative_at_center(pe.prices) - first_derivative_at_center(pe.controlPrices) + delta(black, spot)
  opt.results.gamma = second_derivative_at_center(pe.prices) - second_derivative_at_center(pe.prices) + gamma(black, spot)

  return pe, opt
end
