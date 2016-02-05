using Base.Dates

abstract DateGenerationRule
type DateGenerationBackwards <: DateGenerationRule end
type DateGenerationForwards <: DateGenerationRule end

type Schedule{B <: BusinessDayConvention, D <: DateGenerationRule, C <: BusinessCalendar}
  effectiveDate::Date
  terminationDate::Date
  tenor::TenorPeriod
  convention::B
  termDateConvention::B
  rule::D
  endOfMonth::Bool
  dates::Vector{Date}
  cal::C

  function Schedule(effectiveDate::Date, terminationDate::Date, tenor::TenorPeriod, convention::B, termDateConvention::B,
                  rule::D, endOfMonth::Bool, dates::Vector{Date}, cal::C = TargetCalendar())
    # adjust end date if necessary
    dates[end] = adjust(cal, termDateConvention, dates[end])

    new(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
  end
end

function Schedule{B <: BusinessDayConvention, C <: BusinessCalendar}(effectiveDate::Date, terminationDate::Date, tenor::TenorPeriod, convention::B, termDateConvention::B,
  rule::DateGenerationForwards, endOfMonth::Bool, cal::C = TargetCalendar())
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

  return Schedule{B, DateGenerationForwards, C}(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
end

function Schedule{B <: BusinessDayConvention, C <: BusinessCalendar}(effectiveDate::Date, terminationDate::Date, tenor::TenorPeriod, convention::B, termDateConvention::B,
  rule::DateGenerationBackwards, endOfMonth::Bool, cal::C = TargetCalendar())
  size = get_size(tenor.period, effectiveDate, terminationDate)
  dates = Vector{Date}(size)
  dates[1] = effectiveDate

  dates[end] = terminationDate
  period = 1
  for i = size - 1:-1:2
    dates[i] = adjust(cal, convention, terminationDate - period * tenor.period)
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

  return Schedule{B, DateGenerationBackwards, C}(effectiveDate, terminationDate, tenor, convention, termDateConvention, rule, endOfMonth, dates, cal)
end

# helpers
function get_size(p::Base.Dates.Month, ed::Date, td::Date)
  return Int(ceil(ceil(Int(td - ed) / 30) / Int(p)))
end

function get_size(p::Base.Dates.Year, ed::Date, td::Date)
  # ed_day, ed_month = monthday(ed)
  # td_day, td_month = monthday(td)
  if monthday(ed) == monthday(td)
    return Int(round(Int(td - ed) / 365) + 1)
  else
    return Int(ceil(Int(td - ed) / 365) + 1)
  end
end
