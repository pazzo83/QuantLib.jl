Index Types and Functions
=========================

These types and methods help to construct indexes (e.g. Libor) and use them in instrument pricing.

General Interest Rate Index Methods
-----------------------------------

InterestRateIndex is the Abstract Type from which the Ibor and Libor indexes are derived.  These are some general methods for all interest rate indexes.

.. function:: fixing_date(idx::InterestRateIndex, d::Date)

    Returns the fixing date of the index

.. function:: value_date(idx::InterestRateIndex, d::Date)

    Returns the value date of the index

.. function:: fixing(idx::InterestRateIndex, ts::TermStructure, _fixing_date::Date, forecast_todays_fixing::Bool=true)

    Returns the fixing of the index at a given date

.. function:: forecast_fixing(idx::InterestRateIndex, ts::TermStructure, _fixing_date::Date)

    Calculates a forecasted fixing for a future date, given a fixing date

.. function:: forecast_fixing(idx::InterestRateIndex, ts::TermStructure, d1::Date, d2::Date, t::Float64)

    Calculates a forecasted fixing for a future date given two dates and the time period between the two dates

.. function:: is_valid_fixing_date(idx::InterestRateIndex, d::Date)

    Determines whether or not a date is a valid fixing date for the index (namely, is it a valid business day)

.. function:: add_fixing!(idx::InterestRateIndex, d::Date, fixingVal::Float64)

    Adds a fixing and date to the index's cache of fixings

Ibor Index
----------

.. code-block:: julia

    immutable IborIndex <: InterestRateIndex
      familyName::AbstractString
      tenor::TenorPeriod
      fixingDays::Int
      currency::AbstractCurrency
      fixingCalendar::BusinessCalendar
      convention::BusinessDayConvention
      endOfMonth::Bool
      dc::DayCount
      ts::TermStructure
      pastFixings::Dict{Date, Float64}
    end

.. function:: IborIndex(familyName::AbstractString, tenor::TenorPeriod, fixingDays::Int, currency::AbstractCurrency, fixingCalendar::BusinessCalendar, convention::BusinessDayConvention, endOfMonth::Bool, dc::DayCount, ts::TermStructure = NullTermStructure())

    Constructor for the Ibor Index, will default to a NullTermStructure (which must be set later - this actually will clone this index with a new TS, since the type is immutable)

.. function:: maturity_date(idx::IborIndex, d::Date)

    Returns the maturity date of the index


Libor Index
-----------

.. code-block:: julia

    immutable LiborIndex <: InterestRateIndex
      familyName::AbstractString
      tenor::TenorPeriod
      fixingDays::Int
      currency::Currency
      fixingCalendar::BusinessCalendar
      jointCalendar::JointCalendar
      convention::BusinessDayConvention
      endOfMonth::Bool
      dc::DayCount
      ts::TermStructure
    end

.. function:: LiborIndex(familyName::AbstractString, tenor::TenorPeriod, fixingDays::Int, currency::Currency, fixingCalendar::BusinessCalendar, jointCalendar::JointCalendar, convention::BusinessDayConvention, endOfMonth::Bool, dc::DayCount, ts::TermStructure = NullTermStructure())

    Default constructor for a Libor Index

.. function:: LiborIndex(familyName::AbstractString, tenor::TenorPeriod, fixingDays::Int, currency::Currency, fixingCalendar::BusinessCalendar, dc::DayCount, yts::YieldTermStructure)

    Additional constructor for a Libor Index, with no joint calendar or business day convention passed (the joint calendar is calculated)

.. function:: value_date(idx::LiborIndex, d::Date)

    Returns the value date of a libor index

.. function:: maturity_date(idx::LiborIndex, d::Date)

    Returns the maturity date of a libor index


Indexes Derived from Libor and Ibor
-----------------------------------

.. function:: euribor_index(tenor::TenorPeriod)

    Builds a Euribor Ibor index with a given time period (e.g. 6 month)

.. function:: euribor_index(tenor::TenorPeriod, ts::TermStructure)

    Builds a Euribor Ibor index with a given time period (e.g. 6 month) with a custom term structure

.. function:: usd_libor_index(tenor::TenorPeriod, yts::YieldTermStructure)

    Builds a USD Libor index with a given time period (e.g. 6 month) and term structure


Swap Index
----------

.. code-block:: julia

    immutable SwapIndex <: InterestRateIndex
      familyName::AbstractString
      tenor::TenorPeriod
      fixingDays::Int
      currency::Currency
      fixingCalendar::BusinessCalendar
      fixedLegTenor::Dates.Period
      fixedLegConvention::BusinessDayConvention
      fixedLegDayCount::DayCount
      discount::TermStructure
      iborIndex::IborIndex
      exogenousDiscount::Bool
      lastSwap::VanillaSwap
      lastFixingDate::Date
    end

.. function:: SwapIndex(familyName::AbstractString, tenor::TenorPeriod, fixingDays::Int, currency::Currency, fixingCalendar::BusinessCalendar, fixedLegTenor::Dates.Period, fixedLegConvention::BusinessDayConvention, fixedLegDayCount::DayCount, discount::TermStructure, iborIndex::IborIndex, exogenousDiscount::Bool = true)

    Constructor for a Swap Index

.. function EuriborSwapIsdaFixA(tenor::TenorPeriod, forwardingTS::YieldTermStructure, discTS::YieldTermStructure)

    Builds a Euribor Swap Isda Fix A index
