# main term structures
using Dates

struct NullYieldTermStructure <: YieldTermStructure end

# function calculate!{T <: TermStructure}(ts::T, recalculate::Bool=false)
#   if !ts.calculated || recalculate
#     _calculate!(ts)
#   end
#
#   return ts
# end

struct JumpDate
  ts_quote::Quote
  ts_date::Date
end

struct JumpTime
  ts_quote::Quote
  ts_time::Float64
end

discount(yts::YieldTermStructure, date::Date) = discount(yts, time_from_reference(yts, date))

function discount(yts::YieldTermStructure, time_frac::Float64)
  disc = discount_impl(yts, time_frac)
  if isdefined(yts, :jumpTimes)
    if length(yts.jumpTimes) == 0
      return disc
    end

    jump_effect = 1.0
    for jump in yts.jumpTimes
      if jump.ts_time > 0.0 && jump.ts_time < time
        if jump.ts_quote.value > 0.0 && jump.ts_quote.value <= 1.0
          jump_effect *= jump.ts_quote.value
        end
      end
    end

    return jump_effect * disc
  else
    return disc
  end
end

function zero_rate(yts::YieldTermStructure, date::Date, dc::DayCount, comp::CompoundingType, freq::Frequency = Annual())
  if date == yts.referenceDate
    return implied_rate(1.0 / discount(yts, 0.0001), dc, comp, 0.0001, freq)
  else
    return implied_rate(1.0 / discount(yts, date), dc, comp, reference_date(yts), date, freq)
  end
end

function zero_rate(yts::YieldTermStructure, time_frac::Float64, comp::CompoundingType, freq::Frequency = Annual())
  t = time_frac == 0.0 ? 0.0001 : time_frac
  return implied_rate(1.0 / discount(yts, t), yts.dc, comp, t, freq)
end

function forward_rate(yts::YieldTermStructure, date1::Date, date2::Date, dc::DayCount, comp::CompoundingType, freq::Frequency)
  if date1 == date2
    t1 = max(time_from_reference(yts, date1) - 0.0001 / 2.0, 0.0)
    t2 = t1 + 0.0001
    return implied_rate(discount(yts, t1) / discount(yts, d2), dc, comp, 0.0001, freq)
  elseif date1 < date2
    return implied_rate(discount(yts, date1) / discount(yts, date2), dc, comp, date1, date2, freq)
  else
    error("Forward start date must be before forward end date")
  end
end

forward_rate(yts::YieldTermStructure, date::Date, period::Integer, dc::DayCount, comp::CompoundingType, freq::Frequency) = forward_rate(yts, date, date + Dates.Day(period), dc, comp, freq)

function forward_rate(yts::YieldTermStructure, time1::Float64, time2::Float64, comp::CompoundingType, freq::Frequency)
  if time1 == time2
    t1 = max(time1 - 0.0001 / 2.0, 0.0)
    t2 = t1 + 0.0001
    interval, compound = (t2 - t1, discount(yts, t1) / discount(yts, t2))
  else
    interval, compound = (time2 - time1, discount(yts, time1) / discount(yts, time2))
  end

  return implied_rate(compound, yts.dc, comp, interval, freq)
end

## FlatForwardTermStructure
mutable struct FlatForwardTermStructure{B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency} <: YieldTermStructure
  settlementDays::Int
  referenceDate::Date
  calendar::B
  forward::Quote
  dc::DC
  comp::C
  freq::F
  rate::InterestRate{DC, C, F}
  jumpTimes::Vector{JumpTime}
  jumpDates::Vector{JumpDate}
end

function FlatForwardTermStructure(settlement_days::Int, referenceDate::Date, calendar::B, forward::Quote, dc::DC,
                                  comp::C = ContinuousCompounding(), freq::F = QuantLib.Time.Annual()) where {B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency}
  rate = InterestRate(forward.value, dc, comp, freq)
  FlatForwardTermStructure{B, DC, C, F}(settlement_days, referenceDate, calendar, forward, dc, comp, freq, rate, Vector{JumpTime}(0), Vector{JumpDate}(0))
end

FlatForwardTermStructure(referenceDate::Date, calendar::B, forward::Quote, dc::DC, comp::C = ContinuousCompounding(), freq::F = QuantLib.Time.Annual()) where {B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency} =
                        FlatForwardTermStructure{B, DC, C, F}(0, referenceDate, calendar, forward, dc, comp, freq)

FlatForwardTermStructure(settlementDays::Int, calendar::B, forward::Quote, dc::DC, comp::C = ContinuousCompounding(), freq::F = QuantLib.Time.Annual()) where {B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency} =
                        FlatForwardTermStructure{B, DC, C, F}(settlementDays, Date(0), calendar, forward, dc, comp, freq)

FlatForwardTermStructure(referenceDate::Date, forward::Float64, dc::DC) where {DC <: DayCount} =
                        FlatForwardTermStructure{TargetCalendar, DC, ContinuousCompounding, Annual}(0, referenceDate, TargetCalendar(), Quote(forward), dc, ContinuousCompounding(), Annual())

FlatForwardTermStructure(referenceDate::Date, forward::Float64, dc::DC, compounding::C, freq::F) where {DC <: DayCount, C <: CompoundingType, F <: Frequency} =
                        FlatForwardTermStructure{TargetCalendar, DC, C, F}(0, referenceDate, TargetCalendar(), Quote(forward), dc, compounding, freq)

discount_impl(ffts::FlatForwardTermStructure, time_frac::Float64) = discount_factor(ffts.rate, time_frac)

# discount_impl(ffts::FlatForwardTermStructure, time_frac::Float64) = discount_factor(ffts.rate, time_frac)
