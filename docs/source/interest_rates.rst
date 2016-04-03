Interest Rates
==============

Concrete Interest Rate type

.. code-block:: julia

    type InterestRate
      rate::Float64
      dc::DayCount
      comp::CompoundingType
      freq::Frequency
    end

.. function:: discount_factor(ir::InterestRate, time_frac::Float64)

    Discount factor implied by the rate compounded at time time_frac

.. function:: discount_factor(ir::InterestRate, date1::Date, date2::Date, ref_start::Date = Date(), ref_end::Date = Date())

    Discount factor implied by the rate compounded between two dates

.. function:: compound_factor(ir::InterestRate, time_frac::Float64)

    Compound factor implied by the rate compounded at time time_frac

.. function:: compound_factor(ir::InterestRate, date1::Date, date2::Date, ref_start::Date = Date(), ref_end::Date = Date())

    Compound factor implied by the rate compounded between two dates

.. function:: equivalent_rate(ir::InterestRate, comp::CompoundingType, freq::Frequency, time_frac::Float64)

    Equivalent interest rate for a compounding period time_frac.  The resulting InterestRate shares the same implicit day-counting rule of the original InterestRate instance

.. function:: equivalent_rate(ir::InterestRate, result_dc::DayCount, comp::CompoundingType, freq::Frequency, date1::Date, date2::Date, ref_start::Date = Date(), ref_end::Date = Date())

    Equivalent rate for a compounding period between two dates.  The resulting rate is calculated taking the required day-counting rule into account

.. function:: implied_rate(compound::Float64, dc::DayCount, comp::CompoundingType, time_frac::Float64, freq::Frequency)

    Implied interest rate for a given compound factor at a given time.  The resulting InterestRate has the day-counter provided as input

.. function:: implied_rate(compound::Float64, dc::DayCount, comp::CompoundingType, date1::Date, date2::Date, freq::Frequency, ref_start::Date = Date(), ref_end::Date = Date())

    Implied rate for a given compound factor between two dates.  The resulting rate is calculated taking the required day-counting rule into account


Compounding Types
-----------------

.. code-block:: julia

    abstract CompoundingType

    type ContinuousCompounding <: CompoundingType end # exp(r * t)
    type SimpleCompounding <: CompoundingType end     # (1+r*t)
    type CompoundedCompounding <: CompoundingType end # (1 + r)^t
    type SimpleThenCompounded <: CompoundingType end  # Simple up to the first period then Compounded


Duration
--------

.. code-block:: julia

    abstract Duration
    
    type ModifiedDuration <: Duration end
