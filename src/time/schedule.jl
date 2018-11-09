using Dates

# helper methods
function next_twentieth(d::Date)
  res = Date(year(d), month(d), 20)
  if res < d
    res += Month(1)
  end

  # add logic for certain rules passed in
  m = month(res)
  if m % 3 != 0 # not a main IMM month
    skip_ = 3 - m%3
    res += skip_ * Month(1)
  end

  return res
end

abstract type DateGenerationRule end
struct DateGenerationBackwards <: DateGenerationRule end
struct DateGenerationForwards <: DateGenerationRule end
struct DateGenerationTwentieth <: DateGenerationRule end

struct Schedule{B <: BusinessDayConvention, B1 <: BusinessDayConvention, D <: DateGenerationRule, C <: BusinessCalendar}
  effectiveDate::Date
  terminationDate::Date
  tenor::TenorPeriod
  convention::B
  termDateConvention::B1
  rule::D
  endOfMonth::Bool
  dates::Vector{Date}
  cal::C

  function Schedule{B, B1, D, C}(effectiveDate::Date,
                                terminationDate::Date,
                                tenor::TenorPeriod,
                                convention::B,
                                termDateConvention::B1,
                                rule::D,
                                endOfMonth::Bool,
                                dates::Vector{Date},
                                cal::C = TargetCalendar()) where {B, B1, D, C}
    # adjust end date if necessary
    dates[end] = adjust(cal, termDateConvention, dates[end])

    new{B, B1, D, C}(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
  end
end

function Schedule(effectiveDate::Date,
                  terminationDate::Date,
                  tenor::TenorPeriod,
                  convention::B,
                  termDateConvention::B1,
                  rule::DateGenerationForwards,
                  endOfMonth::Bool,
                  cal::C = TargetCalendar()) where {B <: BusinessDayConvention, B1 <: BusinessDayConvention, C <: BusinessCalendar}
  # dt = effectiveDate
  # num_dates = 1
  #
  # while dt < terminationDate
  #   dt += tenor.period
  #   num_dates += 1
  # end
  #
  # dates = Vector{Date}(num_dates)
  #
  # dates[1] = effectiveDate
  # dates[end] = terminationDate
  #
  # dt = effectiveDate + tenor.period
  # i = 2
  # while dt < terminationDate
  #   dates[i] = dt
  #   dt += tenor.period
  #   i += 1
  # end

  # this way is about 5-10 microseconds faster for semiannual freq over 25 years
  dates = Vector{Date}()
  dt = effectiveDate
  push!(dates, adjust(cal, convention, dt))
  dt += tenor.period
  while dt < terminationDate
    push!(dates, adjust(cal, convention, dt))
    dt += tenor.period
  end

  if dates[end] != terminationDate
    push!(dates, terminationDate)
  end

  return Schedule{B, B1, DateGenerationForwards, C}(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
end

function Schedule(effectiveDate::Date,
                  terminationDate::Date,
                  tenor::TenorPeriod,
                  convention::B,
                  termDateConvention::B1,
                  rule::DateGenerationBackwards,
                  endOfMonth::Bool,
                  cal::C = TargetCalendar()) where {B <: BusinessDayConvention, B1 <: BusinessDayConvention, C <: BusinessCalendar}
  size = get_size(tenor.period, effectiveDate, terminationDate)
  dates = Vector{Date}(undef, size)
  dates[1] = effectiveDate

  dates[end] = terminationDate
  period = 1
  @simd for i = size - 1:-1:2
    @inbounds dates[i] = adjust(cal, convention, terminationDate - period * tenor.period)
    period += 1
  end
  # dt = effectiveDate
  # last_date = terminationDate
  # insert!(dates, 1, terminationDate)
  # period = 1
  # while true
  #   temp = last_date - period * tenor.period
  #   if temp < dt
  #     break
  #   end
  #   insert!(dates, 1, temp)
  #   period += 1
  # end
  #
  # insert!(dates, 1, effectiveDate)

  return Schedule{B, B1, DateGenerationBackwards, C}(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
end

function Schedule(effectiveDate::Date,
                  terminationDate::Date,
                  tenor::TenorPeriod,
                  convention::B,
                  termDateConvention::B1,
                  rule::DateGenerationTwentieth,
                  endOfMonth::Bool,
                  cal::C = TargetCalendar()) where {B <: BusinessDayConvention, B1 <: BusinessDayConvention, C <: BusinessCalendar}

  dates = Vector{Date}()
  dt = effectiveDate
  push!(dates, adjust(cal, convention, dt))
  seed = effectiveDate

  # next 20th
  next20th = next_twentieth(effectiveDate)

  if next20th != effectiveDate
    push!(dates, next20th)
    seed = next20th
  end

  seed += tenor.period
  while seed < terminationDate
    push!(dates, adjust(cal, convention, seed))
    seed += tenor.period
  end

  if dates[end] != adjust(cal, convention, terminationDate)
    push!(dates, next_twentieth(terminationDate))
  else
    push!(dates, terminationDate)
  end

  return Schedule{B, B1, DateGenerationTwentieth, C}(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
end

# helpers
function get_size(p::Dates.Month, ed::Date, td::Date)
  return Int(ceil(ceil(Dates.value(td - ed) / 30) / Dates.value(p)))
end

function get_size(p::Dates.Year, ed::Date, td::Date)
  # ed_day, ed_month = monthday(ed)
  # td_day, td_month = monthday(td)
  if monthday(ed) == monthday(td)
    return Int(round(Dates.value(td - ed) / 365) + 1)
  else
    return Int(ceil(Dates.value(td - ed) / 365) + 1)
  end
end
