# Interest Rates

# Compounding Types
abstract CompoundingType
immutable ContinuousCompounding <: CompoundingType end # exp(r * t)
immutable SimpleCompounding <: CompoundingType end     # (1+r*t)
immutable CompoundedCompounding <: CompoundingType end # (1 + r)^t
immutable SimpleThenCompounded <: CompoundingType end

## Durations ##
immutable ModifiedDuration <: Duration end

type InterestRate{DC <: DayCount, C <: CompoundingType, F <: Frequency}
  rate::Float64
  dc::DC
  comp::C
  freq::F
end

discount_factor(ir::InterestRate, time_frac::Float64) = 1.0 / compound_factor(ir, time_frac)

function discount_factor(ir::InterestRate, date1::Date, date2::Date, ref_start::Date = Date(), ref_end::Date = Date())
  date2 < date1 && error("Date1 $date1 later than date2 $date2")

  return discount_factor(ir, year_fraction(ir.dc, date1, date2, ref_start, ref_end))
end

function compound_factor(ir::InterestRate, time_frac::Float64)
  time_frac < 0.0 && error("negative time not allowed!  Time fraction is $time_frac")

  return _compound_factor(ir.comp, ir.rate, time_frac, ir.freq)
end

# private methods to calculate compound factor based on compounding type
_compound_factor(::SimpleCompounding, rate::Float64, time_frac::Float64, ::Frequency) = 1.0 + rate * time_frac
_compound_factor(::CompoundedCompounding, rate::Float64, time_frac::Float64, freq::Frequency) = (1.0 + rate / QuantLib.Time.value(freq)) ^ (QuantLib.Time.value(freq) * time_frac)
_compound_factor(::ContinuousCompounding, rate::Float64, time_frac::Float64, ::Frequency) = exp(rate * time_frac)
_compound_factor(::SimpleThenCompounded, rate::Float64, time_frac::Float64, freq::Frequency) =
  time_frac <= 1.0 ? _compound_factor(SimpleCompounding(), rate, time_frac, freq) : _compound_factor(CompoundedCompounding(), rate, time_frac, freq)

function compound_factor(ir::InterestRate, date1::Date, date2::Date, ref_start::Date = Date(), ref_end::Date = Date())
  date2 < date1 && error("Date1 $date1 later than date2 $date2")

  return compound_factor(ir, year_fraction(ir.dc, date1, date2, ref_start, ref_end))
end

# equivalent rates
equivalent_rate(ir::InterestRate, comp::CompoundingType, freq::Frequency, time_frac::Float64) = implied_rate(compound_factor(ir, time_frac), ir.dc, comp, time_frac, freq)

function equivalent_rate(ir::InterestRate, result_dc::DayCount, comp::CompoundingType, freq::Frequency, date1::Date, date2::Date, ref_start::Date = Date(), ref_end::Date = Date())
  date1 > date2 && error("Date1 ($date1) later than date2 ($date2)")
  compound = compound_factor(ir, year_fraction(ir.dc, date1, date2, ref_start, ref_end))

  return implied_rate(compound, dc, comp, year_fraction(result_dc, date1, date2, ref_start, ref_end), freq)
end

# implied rates
function implied_rate{DC <: DayCount, C <: CompoundingType, F <: Frequency}(compound::Float64, dc::DC, comp::C, time_frac::Float64, freq::F)
  rate = compound == 1.0 ? 0.0 : _implied_rate(comp, compound, time_frac, freq)

  return InterestRate{DC, C, F}(rate, dc, comp, freq)
end

# methods to calcualte the implied rate based on compounding type
_implied_rate(::SimpleCompounding, compound::Float64, time_frac::Float64, ::Frequency) = (compound - 1.0) / time_frac
_implied_rate(::CompoundedCompounding, compound::Float64, time_frac::Float64, freq::Frequency) = (compound ^ (1.0 / QuantLib.Time.value(freq) * time_frac) - 1.0) * QuantLib.Time.value(freq)
_implied_rate(::ContinuousCompounding, compound::Float64, time_frac::Float64, ::Frequency) = log(compound) / time_frac
_implied_rate(::SimpleThenCompounded, compound::Float64, time_frac::Float64, freq::Frequency) =
  time_frac <= 1.0 ? _implied_rate(SimpleCompounding(), compound, time_frac, freq) : _implied_rate(CompoundedCompounding(), compound, time_frac, freq)

function implied_rate(compound::Float64, dc::DayCount, comp::CompoundingType, date1::Date, date2::Date, freq::Frequency, ref_start::Date = Date(), ref_end::Date = Date())
  date1 > date2 && error("Date1 ($date1) later than date2 ($date2)")

  return implied_rate(compound, dc, comp, year_fraction(dc, date1, date2, ref_start, ref_end), freq)
end

# end
