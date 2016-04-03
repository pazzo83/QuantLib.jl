Time Module
===========

QuantLib.jl's Time sub-module provides calendars, day counters, tenor time periods, and other time related types and methods for pricing.


Misc. Time/Date methods
-----------------------

.. function:: within_next_week(d1::Date, d2::Date)

    Determines whether the second date is within the next week from the first date

.. function:: within_next_week(t1::Float64, t2::Float64)

    Determines whether the second year fraction is within the next week from the first year fraction

.. function:: within_previous_week(d1::Date, d2::Date)

    Determines whether the second date is within the previous week from the first date


Business Calendars
------------------

QuantLib.jl has a number of calendars based on region and asset-type.  Some of this functionality is based off of Ito.jl and BusinessDays.jl

All calendars inherit from the abstract type:

.. code-block:: julia

    abstract BusinessCalendar


General Calendars
~~~~~~~~~~~~~~~~~

A Null calendar exists which has no holidays and is used simply for date movement

.. code-block:: julia

    type NullCalendar <: BusinessCalendar end


Also, a Target Calendar is available which has only basic holidays:
* Saturdays
* Sundays
* New Year's Day (Jan 1)
* Good Friday
* Easter Monday
* Labor Day (May 1)
* Christmas (Dec 25)
* Day of Goodwill (Dec 26)

.. code-block:: julia

    type TargetCalendar <: BusinessCalendar end


A Joint Calendar construction exists that combines two calendars

.. code-block:: julia

    type JointCalendar <: BusinessCalendar
      cal1::B
      cal2::C
    end


Additional calendars are organized by geography and asset type


US Calendars
~~~~~~~~~~~~

.. code-block:: julia

    abstract UnitedStatesCalendar <: WesternCalendar

    type USSettlementCalendar <: UnitedStatesCalendar; end
    type USNYSECalendar <: UnitedStatesCalendar; end
    type USNERCCalendar <: UnitedStatesCalendar; end
    type USGovernmentBondCalendar <: UnitedStatesCalendar; end

**USSettlementCalendar** - General settlement calendar

**USNYSECalendar** - New York Stock Exchange calendar

**USNERCCalendar** - North American Energy Reliability Council calendar

**USGovernmentBondCalendar** - US government bond market


UK Calendars
~~~~~~~~~~~~

.. code-block:: julia

    abstract UnitedKingdomCalendar <: WesternCalendar

    type UKSettlementCalendar <: UnitedKingdomCalendar end
    type UKLSECalendar <: UnitedKingdomCalendar end
    type UKLMECalendar <: UnitedKingdomCalendar end

**UKSettlementCalendar** - UK Settlement calendar

**UKLSECalendar** - London Stock Exchange calendar

**UKLMECalendar** - London Metals Exchange calendar


General Calendar methods
~~~~~~~~~~~~~~~~~~~~~~~~

.. function:: is_business_day(cal::BusinessCalendar, dt::Date)

    Determines whether a given date is a business day, depending on the calendar used

.. function:: easter_date(y::Int)

    Returns the date of Easter Sunday based on year provided

.. function:: is_holiday(::BusinessCalendar, dt::Date)

    Determines whether a given date is a holiday, based on the calendar used

.. function:: advance(days::Day, cal::BusinessCalendar, dt::Date, biz_conv::BusinessDayConvention = Following())

    Advance a number of days from the provided date

.. function:: advance(time_period::Union{Week, Month, Year}, cal::BusinessCalendar, dt::Date, biz_conv::BusinessDayConvention = Following())

    Advance a number of weeks, months, or years from the provided date

.. function:: adjust(cal::BusinessCalendar, d::Date)

    Adjust a date based on a Following business day convention

.. function:: adjust(cal::BusinessCalendar, ::BusinessDayConvention, d::Date)

    Adjust a date based on the provided business day convention


Business Day Conventions
------------------------

These conventions specify the algorithm used to adjust a date in case it is not a valid business day.

.. code-block:: julia

    abstract BusinessDayConvention

    type Unadjusted <: BusinessDayConvention end
    type ModifiedFollowing <: BusinessDayConvention end
    type Following <: BusinessDayConvention end

**Unadjusted** - Do not adjust

**Modified Following** - Choose the first business day after the given holiday unless it belongs to a different month, in which case choose the first business day before the holiday.

**Following** - Choose the first business day after the given holiday.


Day Counters
------------

These types provide methods for determining the length of a time period according to given market convention, both as a number of days and as a year fraction.

Adopted from Ito.jl and InterestRates.jl

All Day Counters inherit from this abstract type:

.. code-block:: julia

    abstract DayCount


**SimpleDayCount** - Simple day counter for reproducing theoretical calculations.

.. code-block:: julia

    type SimpleDayCount <: DayCount end


Day Counter methods
~~~~~~~~~~~~~~~~~~~

.. function:: day_count(c::DayCount, d_start::Date, d_end::Date)

    Returns the number of days between the two dates based off of the day counter method

.. function:: year_fraction(c::DayCount, d_start::Date, d_end::Date)

    Returns the fraction of year between the two dates based off of the day counter method


General Day Counters
~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type Actual360 <:DayCount ; end
    type Actual365 <: DayCount ; end

**Actual360** - Actual / 360 day count convention

**Actual365** - Actual/365 (Fixed) day count convention


30/360 Day Counters
~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    abstract Thirty360 <:DayCount

    type BondThirty360 <: Thirty360; end
    type EuroBondThirty360 <: Thirty360; end
    type ItalianThirty360 <: Thirty360; end

    typealias USAThirty360 BondThirty360
    typealias EuroThirty360 EuroBondThirty360

**USAThirty360** - 30/360 (Bond Basis)

**EuroThirty360** - 30E/360 (Eurobond basis)

**ItalianThirty360** - 30/360 (Italian)


Actual-Actual Day Counters
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    abstract ActualActual <: DayCount

    type ISMAActualActual <: ActualActual; end
    type ISDAActualActual <: ActualActual; end
    type AFBActualActual <: ActualActual; end

    typealias ActualActualBond ISMAActualActual

**ISMAActualActual** - the ISMA and US Treasury convention, also known as "Actual/Actual (Bond)"

**ISDAActualActual** - the ISDA convention, also known as "Actual/Actual (Historical)", "Actual/Actual", "Act/Act", and according to ISDA also "Actual/365", "Act/365", and "A/365"

**AFBActualActual**  - the AFB convention, also known as "Actual/Actual (Euro)"


Frequency
---------

Frequency of events

.. code-block:: julia

    abstract Frequency

    type NoFrequency <: Frequency end
    type Once <: Frequency end
    type Annual <: Frequency end
    type Semiannual <: Frequency end
    type EveryFourthMonth <: Frequency end
    type Quarterly <: Frequency end
    type Bimonthly <: Frequency end
    type Monthly <: Frequency end
    type EveryFourthWeek <: Frequency end
    type Biweekly <: Frequency end
    type Weekly <: Frequency end
    type Daily <: Frequency end
    type OtherFrequency <: Frequency end

.. function:: value(::Frequency)

    Returns the number of times the event will occur in one year (e.g. 1 for Annual, 2 for Semiannual, 3 for EveryFourthMonth, etc)

.. function:: period(::Frequency)

    Returns the underlying time period of the frequency (e.g. 1 Year for Annual, 6 Months for Semiannual, etc)


Schedule
--------

Payment schedule data structure

.. code-block:: julia

    type Schedule
      effectiveDate::Date
      terminationDate::Date
      tenor::TenorPeriod
      convention::BusinessDayConvention
      termDateConvention::BusinessDayConvention
      rule::DateGenerationRule
      endOfMonth::Bool
      dates::Vector{Date}
      cal::BusinessCalendar
    end


Date Generation methods
~~~~~~~~~~~~~~~~~~~~~~~

These conventions specify the rule used to generate dates in a Schedule.

.. code-block:: julia

    abstract DateGenerationRule

    type DateGenerationBackwards <: DateGenerationRule end
    type DateGenerationForwards <: DateGenerationRule end
    type DateGenerationTwentieth <: DateGenerationRule end


**DateGenerationBackwards** - Backward from termination date to effective date.

**DateGenerationForwards** - Forward from effective date to termination date.

**DateGenerationTwentieth** - All dates but the effective date are taken to be the twentieth of their month (used for CDS schedules in emerging markets.)  The termination date is also modified.


Tenor Period
------------

Data structure for a time period with frequency

.. code-block:: julia

    type TenorPeriod
      period::Dates.Period
      freq::Frequency
    end

.. function:: TenorPeriod(f::Frequency)

    Constructor for a TenorPeriod given a Frequency

.. function:: TenorPeriod(p::Dates.Period)

    Constructor for a TenorPeriod given a Julia Date Period


TimeGrid
--------

A Time-Grid data structure

.. code-block:: julia

    type TimeGrid
      times::Vector{Float64}
      dt::Vector{Float64}
      mandatoryTimes::Vector{Float64}
    end

.. function:: TimeGrid(times::Vector{Float64}, steps::Int)

    Time Grid constructor given a vector of times and a number of steps

.. function:: TimeGrid(endTime::Float64, steps::Int)

    Time Grid constructor given an end time and a number of steps

.. function:: closest_time(tg::TimeGrid, t::Float64)

    Returns the time on the grid closest to the given t

.. function:: return_index(tg::TimeGrid, t::Float64)

    Returns the index i such that grid[i] = t
