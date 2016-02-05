# Day Count (adapted from Ito.jl and InterestRates.jl)
using Base.Dates

abstract DayCount

type Actual360 <:DayCount ; end
type Actual365 <: DayCount ; end
abstract Thirty360 <:DayCount

type BondThirty360 <: Thirty360; end
type EuroBondThirty360 <: Thirty360; end
type ItalianThirty360 <: Thirty360; end

typealias USAThirty360 BondThirty360
typealias EuroThirty360 EuroBondThirty360

abstract ActualActual <:DayCount

type ISMAActualActual <: ActualActual; end
type ISDAActualActual <: ActualActual; end
type AFBActualActual <: ActualActual; end

type SimpleDayCount <: DayCount end

# Day Counting
# default day count method
function day_count(c::EuroBondThirty360, d_start::Date, d_end::Date)
  dd1 = day(d_start)
  dd2 = day(d_end)

  mm1 = month(d_start)
  mm2 = month(d_end)

  yy1 = year(d_start)
  yy2 = year(d_end)

  return 360.0 * (yy2 - yy1) + 30.0 * (mm2 - mm1 - 1) + max(0, 30 - dd1) + min(30, dd2)
end

day_count{C <: DayCount}(c::C, d_start::Date, d_end::Date) = Int(d_end - d_start)

# days per year
days_per_year(::Union{Actual360, Thirty360}) = 360.0
days_per_year(::Actual365) = 365.0

# year fractions
# default
year_fraction(c::SimpleDayCount, d_start::Date, d_end::Date) = year_fraction(c, d_start, d_end, Date(), Date())

year_fraction{C <: DayCount}(c::C, d_start::Date, d_end::Date) = day_count(c, d_start, d_end) / days_per_year(c)

# add'l methods
# year_fraction(c::Union{Actual360, Thirty360, Actual365}, d_start::Date, d_end::Date) = year_fraction(c, d_start, d_end, Date(), Date())
year_fraction(c::Union{Actual360, Thirty360, Actual365}, d_start::Date, d_end::Date, ::Date, ::Date) = year_fraction(c, d_start, d_end)

function year_fraction(::SimpleDayCount, d_start::Date, d_end::Date, ::Date, ::Date)
  dm_start = Dates.Day(d_start)
  dm_end = Dates.Day(d_end)

  if dm_start == dm_end || (dm_start > dm_end && Dates.lastdayofmonth(d_end) == d_end) || (dm_start < dm_end && Dates.lastdayofmonth(d_start) == d_start)
    return Int(Dates.Year(d_end) - Dates.Year(d_start)) +  (Int(Dates.Month(d_end)) - Int(Dates.Month(d_start))) / 12.0
  else
    return year_fraction(BondThirty360(), d_start, d_end)
  end
end

function year_fraction(dc::ISDAActualActual, d1::Date, d2::Date, ::Date, ::Date)
  if d1 == d2
    return 0.0
  end

  if d1 > d2
    return -year_fraction(dc, d2, d1, Date(), Date())
  end

  y1 = year(d1)
  y2 = year(d2)

  dib1 = daysinyear(d1)
  dib2 = daysinyear(d2)
  # println(y1)
  # println(y2)
  # println(dib1)
  # println(dib2)

  sum = y2 - y1 - 1

  sum += day_count(dc, d1, Date(y1+1, 1, 1)) / dib1
  # println(d2)
  # println(day_count(dc, Date(y2, 1, 1), d2))
  sum += day_count(dc, Date(y2, 1, 1), d2) / dib2
  return sum
end

function year_fraction(dc::ISMAActualActual, d1::Date, d2::Date, d3::Date, d4::Date)
  if d1 == d2
    return 0.0
  end

  if d1 > d2
    return -year_fraction(d2, d1, d3, d4)
  end

  ref_period_start = d3 != Date() ? d3 : d1
  ref_period_end = d4 != Date() ? d4 : d2

  months = floor(Int, 0.5 + 12 * Int(ref_period_end - ref_period_start) / 365)

  if months == 0
    ref_period_start = d1
    ref_period_end = d1 + Year(1)
    months = 12
  end

  period = months / 12.0

  if d2 <= ref_period_end
    if d1 >= ref_period_start
      return period * day_count(dc, d1, d2) / day_count(dc, ref_period_start, ref_period_end)
    else
      previous_ref = ref_period_start - Month(months)
      if d2 > ref_period_start
        return year_fraction(dc, d1, ref_period_start, previous_ref, ref_period_end) + year_fraction(dc, d1, d2, previous_ref, ref_period_start)
      else
        return year_fraction(dc, d1, d2, previous_ref, ref_period_start)
      end
    end
  else
    sum = year_fraction(dc, d1, ref_period_end, ref_period_start, ref_period_end)
    i = 0
    new_ref_start = new_ref_end = Date()
    while true
      new_ref_start = ref_period_end + Month(i * months)
      new_ref_end = ref_period_end + Month((i + 1) * months)
      if d2 < new_ref_end
        break
      else
        sum += period
        i += 1
      end
    end
    sum += year_fraction(dc, new_ref_start, d2, new_ref_start, new_ref_end)
    return sum
  end
end
