type PiecewiseDefaultCurve{B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap} <: InterpolatedDefaultProbabilityCurve{P}
  lazyMixin::LazyMixin
  settlementDays::Int
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

function perform_calculations!(curve::InterpolatedDefaultProbabilityCurve)
  _calculate!(curve.boot, curve)
  return curve
end

survival_probability(ts::AbstractDefaultProbabilityTermStructure, d::Date) = survival_probability(ts, time_from_reference(ts, d))

function survival_probability(ts::AbstractDefaultProbabilityTermStructure, t::Float64)
  # TODO handle jumps
  return survival_probability_impl(ts, t)
end

function survival_probability_impl(ts::PiecewiseDefaultCurve, t::Float64)
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
default_probability(ts::AbstractDefaultProbabilityTermStructure, t::Float64) = 1.0 - survival_probability(ts, t)

function default_probability(ts::AbstractDefaultProbabilityTermStructure, d1::Date, d2::Date)
  p1 = d1 < reference_date(ts) ? 0.0 : default_probability(ts, d1)
  p2 = default_probability(ts, d2)

  return p2 - p1
end

function default_probability(ts::AbstractDefaultProbabilityTermStructure, t1::Float64, t2::Float64)
  t1 <= t2 || error("time mismatch")

  p1 = t1 < 0.0 ? 0.0 : default_probability(ts, t1)
  p2 = default_probability(ts, t2)
  return p2 - p1
end

default_density(ts::AbstractDefaultProbabilityTermStructure, d::Date) = default_density(ts, time_from_reference(ts, d))
function default_density(ts::AbstractDefaultProbabilityTermStructure, t::Float64)
  # TODO check range
  return default_density_impl(ts, t)
end

default_density_impl(ts::AbstractDefaultProbabilityTermStructure, t::Float64) = hazard_rate_impl(ts, t) * survival_probability_impl(ts, t)
default_density_impl(ts::PiecewiseDefaultCurve, t::Float64) = _default_density_impl(ts, t) * survival_probability_impl(ts, t)

# sub method
function _default_density_impl(ts::PiecewiseDefaultCurve, t::Float64)
  calculate!(ts)
  if t <= ts.times[end]
    return QuantLib.Math.value(ts.interp, t)
  end

  return ts.data[end]
end

hazard_rate(ts::AbstractDefaultProbabilityTermStructure, d::Date) = hazard_rate(ts, time_from_reference(ts, d))

function hazard_rate(ts::AbstractDefaultProbabilityTermStructure, t::Float64)
  S = survival_probability(ts, t)
  return S == 0.0 ? 0.0 : default_density(ts, t) / S
end

function nodes(ts::PiecewiseDefaultCurve)
  calculate!(ts)
  results = Vector{Tuple}(length(ts.dates))
  for i in eachindex(ts.dates)
    results[i] = (ts.dates[i], ts.data[i])
  end

  return results
end
