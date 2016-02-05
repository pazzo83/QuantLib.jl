type TenorPeriod
  period::Dates.Period
  freq::Frequency
end

# Constructors
function TenorPeriod{F <: Frequency}(f::F)
  freq = value(f)
  if freq < 0
    period = Dates.Day(0)
  elseif freq <= 1
    period = Dates.Year(freq)
  elseif freq <= 12
    period = Dates.Month(12 / freq)
  elseif freq <= 52
    period = Dates.Week(52/freq)
  else
    period = Dates.Day(1)
  end

  return TenorPeriod(period, f)
end

TenorPeriod(p::Dates.Year) = TenorPeriod(p, Annual())

function TenorPeriod(p::Dates.Month)
  x = Int(p)
  if x == 6
    return TenorPeriod(p, Semiannual())
  elseif x == 3
    return TenorPeriod(p, Quarterly())
  elseif x == 4
    return TenorPeriod(p, EveryFourthMonth())
  elseif x == 6
    return TenorPeriod(p, Bimonthly())
  else
    return TenorPeriod(p, Monthly())
  end
end

function TenorPeriod(p::Dates.Week)
  x = Int(p)
  if x == 26
    return TenorPeriod(p, Biweekly())
  elseif x == 13
    return TenorPeriod(p, EveryFourthWeek())
  else
    return TenorPeriod(p, Weekly())
  end
end
