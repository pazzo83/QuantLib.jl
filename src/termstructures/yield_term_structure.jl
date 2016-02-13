# main term structures

type NullYieldTermStructure <: YieldTermStructure end

# function calculate!{T <: TermStructure}(ts::T, recalculate::Bool=false)
#   if !ts.calculated || recalculate
#     _calculate!(ts)
#   end
#
#   return ts
# end

type JumpDate
  ts_quote::Quote
  ts_date::Date
end

type JumpTime
  ts_quote::Quote
  ts_time::Float64
end

discount{T <: YieldTermStructure}(yts::T, date::Date) = discount(yts, time_from_reference(yts, date))

function discount{T <: YieldTermStructure}(yts::T, time_frac::Float64)
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

function zero_rate{T <: YieldTermStructure, C <: CompoundingType, F <: Frequency}(yts::T, date::Date, dc::DayCount, comp::C, freq::F)
  if date == yts.referenceDate
    return implied_rate(1.0 / discount(yts, 0.0001), dc, comp, 0.0001, freq)
  else
    return implied_rate(1.0 / discount(yts, date), dc, comp, date, freq)
  end
end

function zero_rate{T <: YieldTermStructure, C <: CompoundingType, F <: Frequency}(yts::T, time_frac::Float64, comp::C, freq::F)
  t = time_frac == 0.0 ? 0.0001 : time_frac
  return implied_rate(1.0 / discount(yts, t), comp, t, freq)
end

function forward_rate{T <: YieldTermStructure, DC <: DayCount, C <: CompoundingType, F <: Frequency}(yts::T, date1::Date, date2::Date, dc::DC, comp::C, freq::F)
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

forward_rate{T <: YieldTermStructure, DC <: DayCount, C <: CompoundingType, F <: Frequency}(yts::T, date::Date, period::Integer, dc::DC, comp::C, freq::F) = forward_rate(yts, date, date + Dates.Day(period), dc, comp, freq)

function forward_rate{T <: YieldTermStructure, C <: CompoundingType, F <: Frequency}(yts::T, time1::Float64, time2::Float64, comp::C, freq::F)
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
type FlatForwardTermStructure{I <: Integer, B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency} <: YieldTermStructure
  settlementDays::I
  referenceDate::Date
  calendar::B
  forward::Quote
  dc::DC
  comp::C
  freq::F
  rate::InterestRate
  jumpTimes::Vector{JumpTime}
  jumpDates::Vector{JumpDate}
end

function FlatForwardTermStructure{I <: Integer, B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency}(settlement_days::I, referenceDate::Date, calendar::B, forward::Quote, dc::DC,
                                  comp::C = ContinuousCompounding(), freq::F = JQuantLib.Time.Annual())
  rate = InterestRate(forward.value, dc, comp, freq)
  FlatForwardTermStructure(settlement_days, referenceDate, calendar, forward, dc, comp, freq, rate, Vector{JumpTime}(0), Vector{JumpDate}(0))
end

FlatForwardTermStructure{B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency}(referenceDate::Date, calendar::B, forward::Quote, dc::DC, comp::C = ContinuousCompounding(), freq::F = JQuantLib.Time.Annual()) =
                        FlatForwardTermStructure(0, referenceDate, calendar, forward, dc, comp, freq)

FlatForwardTermStructure{B <: BusinessCalendar, DC <: DayCount, C <: CompoundingType, F <: Frequency}(settlementDays::Int, calendar::B, forward::Quote, dc::DC, comp::C = ContinuousCompounding(), freq::F = JQuantLib.Time.Annual()) =
                        FlatForwardTermStructure(settlementDays, Date(), calendar, forward, dc, comp, freq)

discount_impl(ffts::FlatForwardTermStructure, time_frac::Float64) = discount_factor(ffts.rate, time_frac)

# discount_impl(ffts::FlatForwardTermStructure, time_frac::Float64) = discount_factor(ffts.rate, time_frac)
