type BlackConstantVol{I <: Integer, BC <: BusinessCalendar, DC <: DayCount} <: BlackVolTermStructure
  referenceDate::Date
  settlementDays::I
  calendar::BC
  volatility::Quote
  dc::DC
end

BlackConstantVol(refDate::Date, cal::BusinessCalendar, volatility::Float64, dc::DayCount) = BlackConstantVol(refDate, 0, cal, Quote(volatility), dc)

black_vol_impl(bts::BlackConstantVol, ::Float64, ::Float64) = bts.volatility.value

function black_variance_impl(bts::BlackConstantVol, t::Float64, strike::Float64)
  vol = black_vol_impl(bts, t, strike)
  return vol * vol * t
end

function black_vol(bts::BlackVolTermStructure, t::Float64, strike::Float64)
  #TODO check stuff
  return black_vol_impl(bts, t, strike)
end

black_vol(bts::BlackVolTermStructure, d::Date, strike::Float64) = black_vol(bts, time_from_reference(bts, d), strike)

function black_variance(bts::BlackVolTermStructure, d::Date, strike::Float64)
  #TODO check stuff
  t = time_from_reference(bts, d)
  return black_variance_impl(bts, t, strike)
end

function black_variance(bts::BlackVolTermStructure, t::Float64, strike::Float64)
  #TODO check stuff
  return black_variance_impl(bts, t, strike)
end
