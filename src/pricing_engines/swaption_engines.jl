## Swaption Pricing Engines ##
using Distributions

type NullSwaptionEngine <: PricingEngine end

_calculate!(::NullSwaptionEngine, ::Swaption) = error("Must set valid pricing engine, use update_pricing_engine")

type BlackSwaptionEngine{Y <: YieldTermStructure, S <: SwaptionVolatilityStructure, DC <: DayCount} <: PricingEngine
  yts::Y
  vol::Quote
  volStructure::S
  dc::DC
  displacement::Float64
end

BlackSwaptionEngine{Y <: YieldTermStructure, DC <: DayCount}(yts::Y, vol::Quote, dc::DC, displacement::Float64 = 0.0) =
                    BlackSwaptionEngine(yts, vol, ConstantSwaptionVolatility(0, JQuantLib.Time.NullCalendar(), JQuantLib.Time.Following(), vol, dc), dc, displacement)

type G2SwaptionEngine{AffineModelType, Y <: YieldTermStructure, I <: Integer} <: PricingEngine
  model::G2{AffineModelType, Y}
  range::Float64
  intervals::I

  call{Y, I}(::Type{G2SwaptionEngine}, model::G2{AffineModelType, Y}, range::Float64, intervals::I) = new{AffineModelType, Y, I}(model, range, intervals)
end

type JamshidianSwaptionEngine{S <: ShortRateModel, Y <: YieldTermStructure} <: PricingEngine
  model::S
  ts::Y

  function call{S}(::Type{JamshidianSwaptionEngine}, model::S)
    new{S, YieldTermStructure}(model)
  end

  function call{S, Y}(::Type{JamshidianSwaptionEngine}, model::S, yts::Y)
    new{S, Y}(model, yts)
  end
end

type LatticeShortRateModelEngineCommon{T <: ShortRateTree}
  tg::TimeGrid
  lattice::T
end


type TreeSwaptionEngine{S <: ShortRateModel, I <: Integer, Y <: YieldTermStructure} <: LatticeShortRateModelEngine{S, Y}
  model::S
  timeSteps::I
  common::LatticeShortRateModelEngineCommon
  ts::Y

  function call{S, I}(::Type{TreeSwaptionEngine}, m::S, tsteps::I)
    t = new{S, I, YieldTermStructure}(m, tsteps)
    add_observer!(m, t)

    return t
  end

  call{S, I}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, l::LatticeShortRateModelEngineCommon) = new{S, I, YieldTermStructure}(m, tsteps, l)

  call{S, I, Y}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, l::LatticeShortRateModelEngineCommon, ts::Y) = new{S, I, T, Y}(m, tsteps, l, ts)
end

function TreeSwaptionEngine{S <: ShortRateModel}(model::S, tg::TimeGrid)
  lattice = tree(model, tg)
  ts = TreeSwaptionEngine(model, 0, LatticeShortRateModelEngineCommon(tg, lattice))

  add_observer!(model, ts)

  return ts
end

type FdG2SwaptionEngine{I <: Integer, Y <: TermStructure} <: PricingEngine
  model::G2
  tGrid::I
  xGrid::I
  yGrid::I
  dampingSteps::I
  invEps::Float64
  schemeDesc::FdmSchemeDesc
  ts::Y
end

FdG2SwaptionEngine(model::G2, tGrid::Int = 100, xGrid::Int = 50, yGrid::Int = 50, dampingSteps::Int = 0, invEps::Float64 = 1e-5,
                  schemeDesc::FdmSchemeDesc = FdmSchemeDesc(Hundsdorfer())) =
                  FdG2SwaptionEngine(model, tGrid, xGrid, yGrid, dampingSteps, invEps, schemeDesc, model.ts)

type FdHullWhiteSwaptionEngine{I <: Integer, Y <: TermStructure} <: PricingEngine
  model::HullWhite
  tGrid::I
  xGrid::I
  dampingSteps::I
  invEps::Float64
  schemeDesc::FdmSchemeDesc
  ts::Y
end

FdHullWhiteSwaptionEngine(model::HullWhite, tGrid::Int = 100, xGrid::Int = 100, dampingSteps::Int = 0, invEps::Float64 = 1e-5,
                  schemeDesc::FdmSchemeDesc = FdmSchemeDesc(Douglas())) =
                  FdHullWhiteSwaptionEngine(model, tGrid, xGrid, dampingSteps, invEps, schemeDesc, model.ts)

# Gaussian One-D engines #
abstract GaussianProbabilities
type NoneProbabilities <: GaussianProbabilities end
type NaiveProbabilities <: GaussianProbabilities end
type DigitalProbabilities <: GaussianProbabilities end

type Gaussian1DSwaptionEngine{G <: Gaussian1DModel, I <: Integer, Y <: YieldTermStructure, P <: GaussianProbabilities} <: PricingEngine
  model::G
  integrationPoints::I
  stddevs::Float64
  extrapolatePayoff::Bool
  flatPayoffExtrapolation::Bool
  discountCurve::Y
  probabilities::P
end

Gaussian1DSwaptionEngine{G <: Gaussian1DModel, I <: Integer, Y <: YieldTermStructure, P <: GaussianProbabilities}(model::G, integrationPoints::I = 64, stddevs::Float64 = 7.0,
                        extrapolatePayoff::Bool = true, flatPayoffExtrapolation::Bool = false, discountCurve::Y = NullYieldTermStructure(), probabilities::P = NoneProbabilities()) =
                        Gaussian1DSwaptionEngine{G, I, Y, P}(model, integrationPoints, stddevs, extrapolatePayoff, flatPayoffExtrapolation, discountCurve, probabilities)

type Gaussian1DNonstandardSwaptionEngine{G <: Gaussian1DModel, I <: Integer, Y <: YieldTermStructure, P <: GaussianProbabilities} <: PricingEngine
  model::G
  integrationPoints::I
  stddevs::Float64
  extrapolatePayoff::Bool
  flatPayoffExtrapolation::Bool
  oas::Quote
  discountCurve::Y
  probabilities::P
end

Gaussian1DNonstandardSwaptionEngine{G <: Gaussian1DModel, I <: Integer, Y <: YieldTermStructure, P <: GaussianProbabilities}(model::G, integrationPoints::I = 64, stddevs::Float64 = 7.0,
                        extrapolatePayoff::Bool = true, flatPayoffExtrapolation::Bool = false, oas::Quote = Quote(-1.0), discountCurve::Y = NullYieldTermStructure(),
                        probabilities::P = NoneProbabilities()) =
                        Gaussian1DNonstandardSwaptionEngine{G, I, Y, P}(model, integrationPoints, stddevs, extrapolatePayoff, flatPayoffExtrapolation, oas, discountCurve, probabilities)

# methods #

function update!(eng::LatticeShortRateModelEngine)
  if length(eng.common.tg.times) > 0
    eng.common.lattice = tree(eng.model, eng.common.tg)
  end

  return eng
end

## General Pricing Functions ##
function black_formula{T <: OptionType}(optionType::T, strike::Float64, forward::Float64, stdDev::Float64, discount::Float64 = 1.0, displacement::Float64 = 0.0)
  # TODO check requirements (see cpp)
  opt_type = JQuantLib.value(optionType)
  if stdDev == 0.0
    return max((forward - strike) * opt_type, 0.0 * discount)
  end

  forward += displacement
  strike += displacement

  if strike == 0.0
    return isa(optionType, Call) ? forward * discount : 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev
  d2 = d1 - stdDev
  norm = Normal() # using distributions.jl
  nd1 = cdf(norm, opt_type * d1)
  nd2 = cdf(norm, opt_type * d2)
  result = discount * opt_type * (forward * nd1 - strike * nd2)

  return result
end

function black_formula_standard_dev_derivative(strike::Float64, forward::Float64, stdDev::Float64, discount::Float64, displacement::Float64)
  forward += displacement
  strike += displacement

  if stdDev == 0.0 || strike == 0.0
    return 0.0
  end

  d1 = log(forward / strike) / stdDev + 0.5 * stdDev

  return discount * forward * distribution_derivative(Normal(), d1)
end

get_annuity(delivery::SettlementPhysical, swap::VanillaSwap) = abs(fixed_leg_BPS(swap)) / basisPoint

function _calculate!(pe::BlackSwaptionEngine, swaption::Swaption)
  exerciseDate = swaption.exercise.dates[1]
  swap = swaption.swap

  strike = swap.fixedRate

  # override swap's pricing engine temporarily, bypassing normal calc flow, since swap.iborIndex might be using a diff curve
  tempSwap = clone(swap, DiscountingSwapEngine(pe.yts))
  # _calculate!(DiscountingSwapEngine(pe.yts), swap)
  calculate!(tempSwap)
  # swap.lazyMixin.calculated = true
  atmForward = fair_rate(tempSwap)

  if tempSwap.spread != 0.0
    correction = tempSwap.spread * abs(floating_leg_BPS(tempSwap) / fixed_leg_BPS(tempSwap))
    strike -= correction
    atmForward -= correction
    swaption.results.additionalResults["spreadCorrection"] = correction
  else
    swaption.results.additionalResults["spreadCorrection"] = 0.0
  end

  swaption.results.additionalResults["strike"] = strike
  swaption.results.additionalResults["atmForward"] = atmForward

  annuity = get_annuity(swaption.delivery, tempSwap)
  swaption.results.additionalResults["annuity"] = annuity

  # the swap length calculation might be improved using the value date of the exercise date
  swapLength = swap_length(pe.volStructure, exerciseDate, date(get_latest_coupon(tempSwap.legs[1])))
  swaption.results.additionalResults["swapLength"] = swapLength

  variance = black_varience(pe.volStructure, exerciseDate, swapLength, strike)

  stdDev = sqrt(variance)
  swaption.results.additionalResults["stdDev"] = stdDev
  w = isa(tempSwap.swapT, Payer) ? Call() : Put()

  swaption.results.value = black_formula(w, strike, atmForward, stdDev, annuity, pe.displacement)

  exerciseTime = time_from_reference(pe.volStructure, exerciseDate)

  swaption.results.additionalResults["vega"] = sqrt(exerciseTime) * black_formula_standard_dev_derivative(strike, atmForward, stdDev, annuity, pe.displacement)

  # # resetting swap
  # reset!(swap.results)
  # swap.lazyMixin.calculated = false

  return swaption
end

function _calculate!(pe::G2SwaptionEngine, swaption::Swaption)
  swap = swaption.swap

  # overriding pricing engine
  _calculate!(DiscountingSwapEngine(pe.model.ts), swap)
  swap.lazyMixin.calculated = true

  correction = swap.spread * abs(floating_leg_BPS(swap) / fixed_leg_BPS(swap))
  fixedRate = swap.fixedRate - correction
  swaption.results.value = gen_swaption(pe.model, swaption, fixedRate, pe.range, pe.intervals)

  return swaption
end

function _calculate!(pe::JamshidianSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model
  ref_date = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  amounts = copy(swaption.swap.args.fixedCoupons)
  amounts[end] += swaption.swap.nominal

  maturity = year_fraction(dc, ref_date, swaption.exercise.dates[1])

  fixedPayTimes = zeros(length(swaption.swap.args.fixedPayDates))
  valueTime = year_fraction(dc, ref_date, swaption.swap.args.fixedResetDates[1])
  for i = 1:length(fixedPayTimes)
    fixedPayTimes[i] = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])
  end

  finder = RStarFinder(tsmodel, swaption.swap.nominal, maturity, valueTime, fixedPayTimes, amounts)

  minStrike = -10.0
  maxStrike = 10.0
  slv = BrentSolver(10000, true, true, minStrike, maxStrike)

  rStar = solve(slv, operator(finder), 1e-8, 0.05, minStrike, maxStrike)

  w = isa(swaption.swap.swapT, Payer) ? Put() : Call()

  _size = length(swaption.swap.args.fixedCoupons)

  val = 0.0
  _B = discount_bond(tsmodel, maturity, valueTime, rStar)

  @simd for i = 1:_size
    @inbounds fixedPayTime = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])

    strike = discount_bond(tsmodel, maturity, fixedPayTime, rStar) / _B

    dboValue = discount_bond_option(tsmodel, w, strike, maturity, valueTime, fixedPayTime)

    @inbounds val += amounts[i] * dboValue
  end

  swaption.results.value = val

  return swaption
end

function greater_than_or_equal_to{T}(x::T, y::T)
  return x >= y
end

function _calculate!(pe::TreeSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model

  refDate = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  dSwaption = DiscretizedSwaption(swaption, refDate, dc)

  if isdefined(pe, :common)
    lattice = pe.common.lattice
  else
    times = mandatory_times(dSwaption)
    tg = TimeGrid(times, pe.timeSteps)
    lattice = tree(pe.model, tg)
    pe.common = LatticeShortRateModelEngineCommon(tg, lattice)
  end

  stoppingTimes = zeros(length(swaption.exercise.dates))
  @simd for i = 1:length(stoppingTimes)
    @inbounds stoppingTimes[i] = year_fraction(dc, refDate, swaption.exercise.dates[i])
  end

  initialize!(dSwaption, lattice.treeLattice, stoppingTimes[end])

  nextExerciseIdx = findnext(greater_than_or_equal_to, stoppingTimes, 1, 0.0)

  nextExercise = nextExerciseIdx != 0 ? stoppingTimes[nextExerciseIdx] : stoppingTimes[end]

  rollback!(dSwaption, nextExercise)
  swaption.results.value = present_value(dSwaption)
end

function _calculate!(pe::FdG2SwaptionEngine, swaption::Swaption)
  # 1. Term structure
  ts = pe.ts

  # 2. Mesher
  dc = ts.dc
  refDate = reference_date(ts)
  maturity = year_fraction(dc, refDate, swaption.exercise.dates[end])

  process1 = OrnsteinUhlenbeckProcess(get_a(pe.model), get_sigma(pe.model))
  process2 = OrnsteinUhlenbeckProcess(get_b(pe.model), get_eta(pe.model))

  xMesher = FdmSimpleProcess1dMesher(pe.xGrid, process1, maturity, 1, pe.invEps)
  yMesher = FdmSimpleProcess1dMesher(pe.yGrid, process2, maturity, 1, pe.invEps)
  mesher = FdmMesherComposite(xMesher, yMesher)

  # 3. Inner Value calculator
  exerciseDates = swaption.exercise.dates
  t2d = Dict{Float64, Date}()
  for i = 1:length(exerciseDates)
    t = year_fraction(dc, refDate, exerciseDates[i])
    t2d[t] = exerciseDates[i]
  end

  disTs = pe.model.ts
  fwdTs = swaption.swap.iborIndex.ts

  fwdModel = G2(fwdTs, get_a(pe.model), get_sigma(pe.model), get_b(pe.model), get_eta(pe.model), get_rho(pe.model))
  calculator = FdmAffineModelSwapInnerValue(pe.model, fwdModel, swaption.swap, t2d, mesher, 1)

  # 4. Step Conditions
  conditions = vanilla_FdmStepConditionComposite(DividendSchedule(), swaption.exercise, mesher, calculator, refDate, dc)

  # 5. Boundary conditions
  boundaries = FdmBoundaryConditionSet()

  # 6. Solver
  solverDesc = FdmSolverDesc(mesher, boundaries, conditions, calculator, maturity, pe.tGrid, pe.dampingSteps)
  solver = FdmG2Solver(pe.model, solverDesc, pe.schemeDesc)
  swaption.results.value = value_at(solver, 0.0, 0.0)
end

function _calculate!(pe::FdHullWhiteSwaptionEngine, swaption::Swaption)
  # 1. Term structure
  ts = pe.ts

  # 2. Mesher
  dc = ts.dc
  refDate = reference_date(ts)
  maturity = year_fraction(dc, refDate, swaption.exercise.dates[end])

  process = OrnsteinUhlenbeckProcess(get_a(pe.model), get_sigma(pe.model))
  shortRateMesher = FdmSimpleProcess1dMesher(pe.xGrid, process, maturity, 1, pe.invEps)
  mesher = FdmMesherComposite(shortRateMesher)

  # 3. Inner Value calculator
  exerciseDates = swaption.exercise.dates
  t2d = Dict{Float64, Date}()
  for i = 1:length(exerciseDates)
    t = year_fraction(dc, refDate, exerciseDates[i])
    t2d[t] = exerciseDates[i]
  end

  disTs = pe.model.ts
  fwdTs = swaption.swap.iborIndex.ts

  # TODO check that day counts and ref dates match btwn the two term structures
  fwdModel = HullWhite(fwdTs, get_a(pe.model), get_sigma(pe.model))
  calculator = FdmAffineModelSwapInnerValue(pe.model, fwdModel, swaption.swap, t2d, mesher, 1)

  # 4. Step Conditions
  conditions = vanilla_FdmStepConditionComposite(DividendSchedule(), swaption.exercise, mesher, calculator, refDate, dc)

  # 5. Boundary conditions
  boundaries = FdmBoundaryConditionSet()

  # 6. Solver
  solverDesc = FdmSolverDesc(mesher, boundaries, conditions, calculator, maturity, pe.tGrid, pe.dampingSteps)

  solver = FdmHullWhiteSolver(pe.model, solverDesc, pe.schemeDesc)
  swaption.results.value = value_at(solver, 0.0)
end
