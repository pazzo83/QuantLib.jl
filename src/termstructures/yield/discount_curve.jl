type InterpolatedDiscountCurve{DC <: DayCount, P <: Interpolation, B <: BusinessCalendar} <: InterpolatedCurve
  settlementDays::Int
  referenceDate::Date
  dc::DC
  interp::P
  cal::B
  dates::Vector{Date}
  times::Vector{Float64}
  data::Vector{Float64}
end

function InterpolatedDiscountCurve(dates::Vector{Date},
                                  discounts::Vector{Float64},
                                  dc::DayCount,
                                  interpolator::Interpolation)
  idc = InterpolatedDiscountCurve(0, dates[1], dc, interpolator, NullCalendar(), dates, zeros(length(dates)), discounts)

  initialize!(idc)

  return idc
end

function initialize!(idc::InterpolatedDiscountCurve)
  length(idc.dates) == length(idc.data) || error("dates / data mismatch")
  idc.times[1] = 0.0
  @simd for i = 2:length(idc.dates)
    @inbounds idc.times[i] = year_fraction(idc.dc, idc.dates[1], idc.dates[i])
  end

  # initialize interpolator
  Math.initialize!(idc.interp, idc.times, idc.data)

  Math.update!(idc.interp)

  return idc
end

function discount_impl(idc::InterpolatedDiscountCurve, t::Float64)
  if t <= idc.times[end]
    return Math.value(idc.interp, t)
  end

  # flat fwd extrapolation
  tMax = idc.times[end]
  dMax = idc.data[end]
  instFwdMax = Math.derivative(idc.interp, tMax) / dMax
  return dMax * exp(-instFwdMax * (t - tMax))
end
