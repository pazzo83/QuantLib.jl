# Business Calendars (adapted from Ito.jl and BusinessDays.jl)
using Base.Dates

abstract BusinessCalendar

abstract WesternCalendar <: BusinessCalendar
abstract OrthodoxCalendar <: BusinessCalendar

# target calendar
type TargetCalendar <: BusinessCalendar end

# for simply moving foward and backward in time
type NullCalendar <: BusinessCalendar end

type JointCalendar{B <: BusinessCalendar, C <: BusinessCalendar} <: BusinessCalendar
  cal1::B
  cal2::C
end

# US Calendars
abstract UnitedStatesCalendar <: WesternCalendar

type USSettlementCalendar <: UnitedStatesCalendar; end
type USNYSECalendar <: UnitedStatesCalendar; end
type USNERCCalendar <: UnitedStatesCalendar; end
type USGovernmentBondCalendar <: UnitedStatesCalendar; end

# UK Calendars
abstract UnitedKingdomCalendar <: WesternCalendar

type UKSettlementCalendar <: UnitedKingdomCalendar end
type UKLSECalendar <: UnitedKingdomCalendar end
type UKLMECalendar <: UnitedKingdomCalendar end

abstract BusinessDayConvention
type Unadjusted <: BusinessDayConvention end
type ModifiedFollowing <: BusinessDayConvention end
type Following <: BusinessDayConvention end

# easter functions
function easter_rata{I <: Integer}(y::I)

  c::Int64
	e::Int64
	p::Int64

   # Algo R only works after 1582
   if y < 1582
        # Are you using this? Send me a postcard!
        error("Year cannot be less than 1582. Provided: $(y).")
   end

	# Century
   c = div( y , 100) + 1

   # Shifted Epact
   e = mod(14 + 11*(mod(y, 19)) - div(3*c, 4) + div(5+8*c, 25), 30)

   # Adjust Epact
   if (e == 0) || ((e == 1) && ( 10 < mod(y, 19) ))
   	e += 1
   end

   # Paschal Moon
   p = Date(y, 4, 19).instant.periods.value - e

   # Easter: locate the Sunday after the Paschal Moon
   return p + 7 - mod(p, 7)
end

# Returns Date
function easter_date{I <: Integer}(y::I)
	# Compute the gregorian date for Rata Die number
     return Date(Dates.rata2datetime( easter_rata(y) ))
end

# calendar functions
# advance{B <: BusinessDayConvention}(time_period::Day, cal::NullCalendar, dt::Date, ::B) = dt += time_period
# advance{B <: BusinessDayConvention}(time_period::Union{Week, Month, Year}, cal::NullCalendar, dt::Date, ::B) = dt += time_period

function advance{C <: BusinessCalendar, B <: BusinessDayConvention}(days::Day, cal::C, dt::Date, biz_conv::B = Following())
  n = Int(days)
  if n > 0
    while n > 0
      dt += Day(1)
      while !is_business_day(cal, dt)
        dt += Day(1)
      end
      n -= 1
    end
  else
    while (n < 0)
      dt -= Day(1)
      while !is_business_day(cal, dt)
        dt -= Day(1)
      end
      n += 1
    end
  end

  return dt
end

function advance{C <: BusinessCalendar, B <: BusinessDayConvention}(time_period::Union{Week, Month, Year}, cal::C, dt::Date, biz_conv::B = Following())
  dt += time_period
  return adjust(cal, biz_conv, dt)
end

is_business_day(cal::NullCalendar, ::Date) = true

function is_business_day{C <: BusinessCalendar}(cal::C, dt::Date)
  if dayofweek(dt) in [6, 7] || is_holiday(cal, dt)
    return false
  else
    return true
  end
end

# In the United States, if a holiday falls on Saturday, it's observed on the preceding Friday.
# If it falls on Sunday, it's observed on the next Monday.
function adjustweekendholidayUS(dt::Date)
	if dayofweek(dt) == 6
		return dt - Dates.Day(1)
	end

	if dayofweek(dt) == 7
		return dt + Dates.Day(1)
	end

	return dt
end

is_holiday(joint::JointCalendar, dt::Date) = is_holiday(joint.cal1, dt) || is_holiday(joint.cal2, dt)

function is_holiday(::USSettlementCalendar , dt::Date)

	const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)

	if (
			# New Year's Day
			adjustweekendholidayUS(Date(yy, 1, 1)) == dt
			||
			# New Year's Day on the previous year when 1st Jan is Saturday
			(mm == 12 &&  dd == 31 && dayofweek(dt) == Friday)
			||
			# Birthday of Martin Luther King, Jr.
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 1)
			||
			# Washington's Birthday
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 2)
			||
			# Memorial Day is the last Monday in May
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 5)
			||
			# Independence Day
			adjustweekendholidayUS(Date(yy, 7, 4)) == dt
			||
			# Labor Day is the first Monday in September
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 1 && mm == 9)
			||
			# Columbus Day is the second Monday in October
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 2 && mm == 10)
			||
			# Veterans Day
			adjustweekendholidayUS(Date(yy, 11, 11)) == dt
			||
			# Thanksgiving Day is the fourth Thursday in November
			(dayofweek(dt) == 4 && dayofweekofmonth(dt) == 4 && mm == 11)
			||
			# Christmas
			adjustweekendholidayUS(Date(yy, 12, 25)) == dt
		)
		return true
	end

	return false
end

function is_holiday(::USGovernmentBondCalendar, dt::Date)
  const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)
	if (
			# New Year's Day
			adjustweekendholidayUS(Date(yy, 1, 1)) == dt
			||
			# New Year's Day on the previous year when 1st Jan is Saturday
			(mm == 12 &&  dd == 31 && dayofweek(dt) == Friday)
			||
			# Birthday of Martin Luther King, Jr.
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 1)
			||
			# Washington's Birthday
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 2)
			||
      # Good Friday
      easter_date(yy) - Day(2) == dt
      ||
			# Memorial Day is the last Monday in May
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 5)
			||
			# Independence Day
			adjustweekendholidayUS(Date(yy, 7, 4)) == dt
			||
			# Labor Day is the first Monday in September
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 1 && mm == 9)
			||
			# Columbus Day is the second Monday in October
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 2 && mm == 10)
			||
			# Veterans Day
			adjustweekendholidayUS(Date(yy, 11, 11)) == dt
			||
			# Thanksgiving Day is the fourth Thursday in November
			(dayofweek(dt) == 4 && dayofweekofmonth(dt) == 4 && mm == 11)
			||
			# Christmas
			adjustweekendholidayUS(Date(yy, 12, 25)) == dt
		)
		return true
	end

	return false
end

function is_holiday(::USNYSECalendar, dt::Date)
  const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)
	if (
			# New Year's Day
			adjustweekendholidayUS(Date(yy, 1, 1)) == dt
			||
			# New Year's Day on the previous year when 1st Jan is Saturday
			(mm == 12 &&  dd == 31 && dayofweek(dt) == Friday)
			||
      # Birthday of Martin Luther King, Jr.
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 1)
			||
			# Washington's Birthday
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) ==3 && mm == 2)
			||
      # Good Friday
      easter_date(yy) - Day(2) == dt
      ||
			# Memorial Day is the last Monday in May
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 5)
			||
			# Independence Day
			adjustweekendholidayUS(Date(yy, 7, 4)) == dt
			||
			# Labor Day is the first Monday in September
			(dayofweek(dt) == 1 && dayofweekofmonth(dt) == 1 && mm == 9)
			||
			# Thanksgiving Day is the fourth Thursday in November
			(dayofweek(dt) == 4 && dayofweekofmonth(dt) == 4 && mm == 11)
			||
			# Christmas
			adjustweekendholidayUS(Date(yy, 12, 25)) == dt
		)
		return true
	end

  # Special Closings
  if (
    # Hurricane Sandy
    (yy == 2012 && mm == 10 && dd in (29, 30))
    ||
    # President Ford's funeral
    (yy == 2007 && mm == 1 && dd == 2)
    ||
    # President Reagan's funeral
    (yy == 2004 && mm == 6 && dd == 11)
    ||
    # 9/11
    (yy == 2001 && mm == 9 && dd in (11, 12, 13, 14))
    )
    return true
  end

	return false
end

function adjustweekendholidayUK(dt::Date)

  if dayofweek(dt) == 6
		return dt + Day(2)
	end

  if dayofweek(dt) == 7
		return dt + Day(1)
	end

	return dt
end

## UK Calendar functions
function is_holiday(c::Union{UKSettlementCalendar, UKLSECalendar}, dt::Date)
  const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)

	if (
    # New Year's Day
    adjustweekendholidayUK(Date(yy, 1, 1)) == dt
    ||
    # Good Friday
    easter_date(yy) - Day(2) == dt
    ||
    # Easter MONDAY
    easter_date(yy) + Day(1) == dt
    ||
    # first MONDAY of May (Early May Bank Holiday)
    (dayofweek(dt) == 1 && dayofweekofmonth(dt) == 1 && mm == 5)
    ||
    # last MONDAY of MAY (Spring Bank Holiday)
    (dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 5)
    ||
    # last MONDAY of August (Summer Bank Holiday)
    (dayofweek(dt) == 1 && dayofweekofmonth(dt) == daysofweekinmonth(dt) && mm == 8)
    ||
    # Christmas (possibly moved to MONDAY or Tuesday)
    adjustweekendholidayUK(Date(yy, 12, 25)) == dt
    ||
    # Boxing Day (possibly moved to MONDAY or TUESDAY)
    adjustweekendholidayUK(adjustweekendholidayUK(Date(yy, 12, 25)) + Day(1)) == dt
    )
    return true
  end

  # Fixed holidays
  if (
    # Substitute date for Spring Bank Holiday
		(dt == Date(2012, 06, 04))
		||
		# Diamond Jubilee of Queen Elizabeth II.
		(dt == Date(2012, 06, 05))
		||
		# Golden Jubilee of Queen Elizabeth II.
		(dt == Date(2002, 06, 03))
		||
		# Substitute date for Spring Bank Holiday
		(dt == Date(2002, 06, 04))
		||
		# Wedding of Prince William and Catherine Middleton
		(dt == Date(2011, 04, 29))
    )
    return true
  end

  return false
end

function is_holiday(::TargetCalendar, dt::Date)
  const yy = year(dt)
	const mm = month(dt)
	const dd = day(dt)
  const easter_sun = easter_date(yy)

  if (
    # New Years Day
    (mm == 1 && dd == 1)
    ||
    # Good Friday
    easter_sun - Day(2) == dt
    ||
    # Easter Monday
    easter_sun + Day(1) == dt
    ||
    # Int'l Labor Day
    (mm == 5 && dd == 1)
    ||
    # Christmas
    (mm == 12 && dd == 25)
    ||
    # Day of Goodwill
    (mm == 12 && dd == 26)
    )
    return true
  else
    return false
  end
end

# adjustments
adjust{B <: BusinessCalendar}(::B, ::Unadjusted, d::Date) = d

function adjust{B <: BusinessCalendar}(cal::B, ::Union{ModifiedFollowing, Following}, d::Date)
  while !is_business_day(cal, d)
    d += Day(1)
  end

  return d
end
