type InterpolatedHazardRateCurve{DC <: DayCount, P <: Interpolation, B <: BusinessCalendar} <: InterpolatedDefaultProbabilityCurve{P}
  settlementDays::Int
  referenceDate::Date
  dc::DC
  interp::P
  cal::B
  dates::Vector{Date}
  times::Vector{Float64}
  data::Vector{Float64}
end

function InterpolatedHazardRateCurve(dates::Vector{Date},
                                    hazardRates::Vector{Float64},
                                    dc::DayCount,
                                    interpolator::Interpolation)
  ihc = InterpolatedHazardRateCurve(0, dates[1], dc, interpolator, NullCalendar(), dates, zeros(length(dates)), hazardRates)

  initialize!(ihc)

  return ihc
end

function initialize!(ihc::InterpolatedHazardRateCurve)
  length(ihc.dates) == length(ihc.data) || error("dates / data mismatch")

  ihc.times[1] = 0.0
  @simd for i = 2:length(ihc.dates)
    @inbounds ihc.times[i] = year_fraction(ihc.dc, ihc.dates[1], ihc.dates[i])
  end

  # initialize interpolator
  Math.initialize!(ihc.interp, ihc.times, ihc.data)

  Math.update!(ihc.interp)

  return ihc
end

function hazard_rate_impl(ihc::InterpolatedHazardRateCurve, t::Float64)
  if t <= ihc.times[end]
    return Math.value(ihc.interp, t)
  end

  # flat hazard rate extrap
  return ihc.data[end]
end

function survival_probability_impl(ihc::InterpolatedHazardRateCurve, t::Float64)
  if t == 0.0
    return 1.0
  end

  if t <= ihc.times[end]
    integral = get_primitive(ihc.interp, t, true)
  else
    # flat hazard rate extrapolation
    integral = get_primitive(ihc.interp, ihc.times[end], true) + ihc.data[end] * (t - ihc.times[end])
  end

  return exp(-integral)
end
