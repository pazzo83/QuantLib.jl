# QuantLib Time module
module Time

# Date helper funcs
within_next_week(d1::Date, d2::Date) = d2 >= d1 && d2 <= d1 + Dates.Day(7)
within_next_week(t1::Float64, t2::Float64) = t1 <= t2 && t2 <= t1 + (1.0/52.0)
within_previous_week(d1::Date, d2::Date) = d2 >= d1 - Dates.Day(7) && d2 <= d1

export within_next_week, within_previous_week

# Frequency.jl
export Frequency, NoFrequency, Once, Annual, Semiannual, EveryFourthMonth, Quarterly, Bimonthly, Monthly, EveryFourthWeek, Biweekly, Weekly, Daily, OtherFrequency,
value, period

# DayCount.jl
export DayCount, Actual360, Actual365, Thirty360, BondThirty360, EuroBondThirty360, ItalianThirty360, ActualActual, SimpleDayCount, ISMAActualActual, ISDAActualActual,
day_count, days_per_year, year_fraction

# BusinessCalendar.jl
export BusinessCalendar, WesternCalendar, OrthodoxCalendar, UnitedStatesCalendar, USSettlementCalendar, USNYSECalendar, USNERCCalendar,
USGovernmentBondCalendar, UnitedKingdomCalendar, UKSettlementCalendar, UKLSECalendar, UKLMECalendar, TargetCalendar, NullCalendar, JointCalendar,
BusinessDayConvention, Unadjusted, ModifiedFollowing, Following,
easter_date, is_holiday, advance, adjust

# tenor_period.jl
export TenorPeriod

# schedule.jl
export DateGenerationRule, DateGenerationForwards, DateGenerationBackwards, DateGenerationTwentieth, Schedule

# time_grid.jl
export TimeGrid, is_empty, closest_time, return_index


include("Frequency.jl")
include("DayCount.jl")
include("BusinessCalendar.jl")
include("tenor_period.jl")
include("schedule.jl")
include("time_grid.jl")

end
