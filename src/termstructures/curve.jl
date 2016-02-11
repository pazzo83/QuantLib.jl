# Curves

type NullCurve <: Curve end

type PiecewiseYieldCurve{I <: Integer, B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap} <: InterpolatedCurve{P, T}
  lazyMixin::LazyMixin
  settlementDays::I
  referenceDate::Date
  instruments::Vector{B}
  dc::DC
  interp::P
  trait::T
  accuracy::Float64
  boot::BT
  times::Vector{Float64}
  dates::Vector{Date}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseYieldCurve{B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap}(referenceDate::Date, instruments::Vector{B}, dc::DC, interp::P, trait::T,
                                accuracy::Float64, boot::BT = IterativeBootstrap())
  # get the initial length of instruments
  n = length(instruments)
  # create an initial state of the curve
  pyc = PiecewiseYieldCurve(LazyMixin(),
                            0,
                            referenceDate,
                            instruments,
                            dc,
                            interp,
                            trait,
                            accuracy,
                            boot,
                            Vector{Float64}(n + 1),
                            Vector{Date}(n + 1),
                            Vector{Float64}(n + 1),
                            Vector{Function}(n + 1),
                            false)

  # initialize the bootstrapping
  initialize(pyc.boot, pyc)

  return pyc
end

type PiecewiseDefaultCurve{I <: Integer, B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap} <: InterpolatedDefaultProbabilityCurve{P, T}
  lazyMixin::LazyMixin
  settlementDays::I
  referenceDate::Date
  instruments::Vector{B}
  dc::DC
  interp::P
  trait::T
  accuracy::Float64
  boot::BT
  times::Vector{Float64}
  dates::Vector{Date}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseDefaultCurve{B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap}(referenceDate::Date, instruments::Vector{B}, dc::DC, interp::P, trait::T,
                                accuracy::Float64, boot::BT = IterativeBootstrap())
  # get the initial length of instruments
  n = length(instruments)
  # create an initial state of the curve
  pyc = PiecewiseDefaultCurve(LazyMixin(),
                            0,
                            referenceDate,
                            instruments,
                            dc,
                            interp,
                            trait,
                            accuracy,
                            boot,
                            Vector{Float64}(n + 1),
                            Vector{Date}(n + 1),
                            Vector{Float64}(n + 1),
                            Vector{Function}(n + 1),
                            false)

  # initialize the bootstrapping
  initialize(pyc.boot, pyc)

  return pyc
end

type FittedBondDiscountCurve{I <: Integer, C <: BusinessCalendar, B <: BondHelper, DC <: DayCount, F <: FittingMethod} <: Curve
  lazyMixin::LazyMixin
  settlementDays::I
  referenceDate::Date
  calendar::C
  bonds::Vector{B}
  dc::DC
  fittingMethod::F
  accuracy::Float64
  maxEvaluations::I
  simplexLambda::Float64

  FittedBondDiscountCurve(settlementDays::I,
                          referenceDate::Date,
                          calendar::C,
                          bonds::Vector{B},
                          dc::DC,
                          fittingMethod::F,
                          accuracy::Float64,
                          maxEvaluations::I,
                          simplexLambda::Float64) =

                          (x = new(LazyMixin(), settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda);
                          x.fittingMethod.commons.costFunction.curve = x)

  # FittedBondDiscountCurve(settlementDays::Integer, referenceDate::Date, calendar::BusinessCalendar, bonds::Vector{B}, dc::DayCount, fittingMethod::FittingMethod, accuracy::Float64=1e-10,
  #                                     maxEvaluations::Integer=10000, simplexLambda::Float64=1.0) =
  #                         new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)
    #n = length(bonds)
  #   println("hi")
  #   this = new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)
  #   println("also hi")
  #   this.fittingMethod.costFunction.curve = this
  #
  #   return this
  # end
end

FittedBondDiscountCurve{I <: Integer, C <: BusinessCalendar, B <: BondHelper, DC <: DayCount, F <: FittingMethod}(settlementDays::I, referenceDate::Date, calendar::C, bonds::Vector{B}, dc::DC, fittingMethod::F, accuracy::Float64=1e-10,
                                     maxEvaluations::I=10000, simplexLambda::Float64=1.0) =
                                     FittedBondDiscountCurve{I, C, B, DC, F}(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)

type FittingCost{I <: Integer} <: CostFunction
  # value::Vector{Float64}
  # values::Vector{Float64}
  firstCashFlow::Vector{I}
  curve::Curve
end

function FittingCost{I <: Integer, C <: Curve}(size::I, curve::C)
  # value = Vector{Float64}()
  # values = Vector{Float64}()
  firstCashFlow = zeros(I, size)

  return FittingCost(firstCashFlow, curve)
end

# Interpolated Curve methods #
max_date{C <: InterpolatedCurve}(curve::C) = curve.dates[end]

function discount{C <: Curve}(curve::C, t::Float64)
  return discount_impl(curve, t)
end

function discount_impl{C <: InterpolatedCurve}(curve::C, t::Float64)
  if t <= curve.times[end]
    return JQuantLib.Math.value(curve.interp, t)
  end

  # println("outside")
  # println(t)
  # println(curve.times)

  # do flat fwd extrapolation
end

function perform_calculations!{C <: InterpolatedCurve}(curve::C)
  _calculate!(curve.boot, curve)
  return curve
end

function perform_calculations!(curve::InterpolatedDefaultProbabilityCurve)
  _calculate!(curve.boot, curve)
  return curve
end


### Fitted curve methods ###
discount_impl(curve::FittedBondDiscountCurve, t::Float64) = discount_function(curve.fittingMethod, curve.fittingMethod.commons.solution, t)

function initialize!(curve::FittedBondDiscountCurve)
  # yield conventions
  dc = curve.dc
  yield_comp = CompoundedCompounding()
  freq = Annual()

  n = length(curve.bonds)
  cost_f = curve.fittingMethod.commons.costFunction

  squared_sum = 0.0
  for i = 1:n
    bond = curve.bonds[i].bond
    leg = bond.cashflows
    clean_price = bond.faceAmount
    bond_settlement = get_settlement_date(bond)

    # get the ytm of the bond
    ytm = yield(bond, clean_price, dc, yield_comp, freq, bond_settlement)
    dur = duration(bond, ytm, dc, yield_comp, freq, ModifiedDuration(), bond_settlement)

    curve.fittingMethod.commons.weights[i] = 1.0 / dur
    squared_sum += curve.fittingMethod.commons.weights[i] * curve.fittingMethod.commons.weights[i]

    cf = bond.cashflows
    for k = 1:length(cf.coupons) # for redemption
      # cf_to_use = k > length(cf.coupons) ? cf.redemption : cf.coupons[i]
      if !has_occurred(cf.coupons[i], bond_settlement)
        cost_f.firstCashFlow[i] = k
        break
      end
    end
  end

  divide_array_by_self!(curve.fittingMethod.commons.weights, sqrt(squared_sum))

  return curve
end

function perform_calculations!(curve::FittedBondDiscountCurve)
  cost_f = curve.fittingMethod.commons.costFunction
  constraint = NoConstraint()

  x = zeros(curve.fittingMethod.size)
  # x = Vector{DD}(curve.fittingMethod.size)

  if length(curve.fittingMethod.commons.guessSolution) > 0
    x = curve.fittingMethod.commons.guessSolution
  end

  simplex = Simplex(curve.simplexLambda)
  problem = Problem(cost_f, constraint, x)

  max_stationary_state_iterations = 100
  root_epsilon = curve.accuracy
  function_epsilon = curve.accuracy
  gradient_norm_epsilon = curve.accuracy

  end_criteria = EndCriteria(curve.maxEvaluations, max_stationary_state_iterations, root_epsilon, function_epsilon, gradient_norm_epsilon)

  minimize!(simplex, problem, end_criteria)
  curve.fittingMethod.commons.solution = problem.currentValue

  number_of_iterations = problem.functionEvaluation
  cost_value = problem.functionValue

  curve.fittingMethod.commons.guessSolution = curve.fittingMethod.commons.solution
  curve.fittingMethod.commons.numberOfIterations = number_of_iterations
  curve.fittingMethod.commons.minimumCostValue = cost_value

  # we have calculated

  return curve
end

function value{C <: CostFunction, T}(cf::C, x::Vector{T})
  ref_date = cf.curve.referenceDate
  dc = cf.curve.dc
  squared_error = 0.0
  n = length(cf.curve.bonds)

  # for (i, bh) in enumerate(cf.curve.bonds)
  for i=1:length(cf.curve.bonds)
    bond = cf.curve.bonds[i].bond
    bond_settlement = get_settlement_date(bond)
    model_price = -accrued_amount(bond, bond_settlement)
    leg = bond.cashflows
    for k = cf.firstCashFlow[i]:length(leg.coupons)
      # @inbounds df = discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.coupons[k])))
      @inbounds model_price += amount(leg.coupons[k]) * discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.coupons[k])))
    end

    # redemption
    # @inbounds model_price += amount(leg.redemption) * discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.redemption)))

    # adjust NPV for forward settlement
    if bond_settlement != ref_date
      model_price /= discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, bond_settlement))
    end

    market_price = bond.faceAmount
    price_error = model_price - market_price
    @inbounds weighted_error = cf.curve.fittingMethod.commons.weights[i] * price_error
    squared_error += weighted_error * weighted_error
  end

  return squared_error
end

survival_probability(ts::AbstractDefaultProbabilityTermStructure, d::Date) = survival_probability(ts, time_from_reference(ts, d))

function survival_probability(ts::AbstractDefaultProbabilityTermStructure, t::Float64)
  # TODO handle jumps
  return survival_probability_impl(ts, t)
end

function survival_probability_impl(ts::InterpolatedHazardRateCurve, t::Float64)
  calculate!(ts)
  if t == 0.0
    return 1.0
  end

  if t <= ts.times[end]
    integral = get_primitive(ts.interp, t, true)
  else
    # flat hazard rate extrapolation
    integral = get_primitive(ts.interp, ts.times[end], true) + ts.data[end] * (t - ts.times[end])
  end

  return exp(-integral)
end

default_probability(ts::AbstractDefaultProbabilityTermStructure, d::Date) = 1.0 - survival_probability(ts, d)

function default_probability(ts::AbstractDefaultProbabilityTermStructure, d1::Date, d2::Date)
  p1 = d1 < reference_date(ts) ? 0.0 : default_probability(ts, d1)
  p2 = default_probability(ts, d2)

  return p2 - p1
end

default_density(ts::AbstractDefaultProbabilityTermStructure, d::Date) = default_density(ts, time_from_reference(ts, d))
function default_density(ts::AbstractDefaultProbabilityTermStructure, t::Float64)
  # TODO check range
  return default_density_impl(ts, t)
end

default_density_impl(ts::InterpolatedHazardRateCurve, t::Float64) = _default_density_impl(ts, t) * survival_probability_impl(ts, t)

# sub method
function _default_density_impl(ts::InterpolatedHazardRateCurve, t::Float64)
  calculate!(ts)
  if t <= ts.times[end]
    return JQuantLib.Math.value(ts.interp, t)
  end

  return ts.data[end]
end

hazard_rate(ts::AbstractDefaultProbabilityTermStructure, d::Date) = hazard_rate(ts, time_from_reference(ts, d))

function hazard_rate(ts::AbstractDefaultProbabilityTermStructure, t::Float64)
  S = survival_probability(ts, t)
  return S == 0.0 ? 0.0 : default_density(ts, t) / S
end

function nodes(ts::InterpolatedHazardRateCurve)
  calculate!(ts)
  results = Vector{Tuple}(length(ts.dates))
  for i in eachindex(ts.dates)
    results[i] = (ts.dates[i], ts.data[i])
  end

  return results
end
