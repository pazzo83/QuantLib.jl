CashFlow Types and Functions
============================

Types and methods to build and work with cashflows.

Basic CashFlow and Coupon Types and Methods
-------------------------------------------
These are common methods for all coupons in QuantLib.jl.  Coupon is an abstract type from which all other coupons are derived.  A coupon itself is a type of cash flow.

Accessor methods
~~~~~~~~~~~~~~~~

.. function:: accrual_start_date(c::Coupon)

    Returns the accrual start date of the coupon

.. function:: accrual_end_date(c::Coupon)

    Returns the accrual end date of the coupon

.. function:: date_accrual_end(c::Coupon)

    Returns the accrual end date of the coupon

.. function:: ref_period_start(c::Coupon)

    Returns the reference period start date of the coupon

.. function:: ref_period_end(c::Coupon)

    Returns the reference period end date of the coupon

.. function:: get_dc(c::Coupon)

    Returns the day counter of the coupon

.. function:: accrual_period!(c::Coupon, val::Float64)

    Sets the accrual period for a coupon

.. function:: accrual_period(c::Coupon)

    Returns the accrual period for the coupon

.. function:: date(c::Coupon)

    Returns the coupon's payment date


Simple Cash Flow
~~~~~~~~~~~~~~~~

This is a basic cash flow type, typically used for redemptions

.. code-block:: julia

    type SimpleCashFlow <: CashFlow
      amount::Float64
      date::Date
    end

.. function:: amount(cf::SimpleCashFlow)

    Returns the simple cash flow amount

.. function:: date(cf::SimpleCashFlow)

    Returns the simple cash flow date

.. function:: date_accrual_end(cf::SimpleCashFlow)

    Returns the accrual end date of the simple cash flow (which is the date)


Dividend
~~~~~~~~

This type is for any dividend cash flow

.. code-block:: julia

    type Dividend <: CashFlow
      amount::Float64
      date::Date
    end

.. function:: amount(cf::Dividend)

    Returns the dividend amount

.. function:: date(cf::Dividend)

    Returns the dividend date

.. function:: date_accrual_end(cf::Dividend)

    Returns the accrual end date of the dividend (which is the date)


Cash Flow Legs
--------------

Cash Flow legs are basically holders of cash flows and coupons

General CashFlow Leg methods
~~~~~~~~~~~~~~~~~~~

.. function:: npv(leg::Leg, yts::YieldTermStructure, settlement_date::Date, npv_date::Date)

    Calculates the npv of a cash flow leg given a term structure, settlement date, and npv date

.. function:: npv(leg::Leg, y::InterestRate, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date)

    Calculates the npv of a cash flow leg given an interest rate, settlement date, and npv date

.. function:: npvbps(leg::Leg, yts::YieldTermStructure, settlement_date::Date, npv_date::Date, includeSettlementDateFlows::Bool = true)

    Calculates the npv and bps of a cash flow leg given a term structure, settlement date, and npv date

.. function:: duration(::ModifiedDuration, leg::Leg, y::InterestRate, dc::DayCount, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date = Date())

    Calculates the modified duration of a cash flow leg given an interest rate, day count, settlement date, and npv date

.. function:: previous_cashflow_date(cf::Leg, settlement_date::Date)

    Returns the previous cash flow date of a cash flow leg given a settlement date

.. function:: accrual_days(cf::CashFlows, dc::DayCount, settlement_date::Date)

    Returns the number of accrued days of a cashflow given a day counter and a settlement date

.. function:: next_cashflow(cf::Leg, settlement_date::Date)

    Returns the next cash flow from a cash flow leg given a settlement date

.. function:: accrued_amount(cf::Leg, settlement_date::Date, include_settlement_cf::Bool = false)

    Calculates the accrued amount from a cash flow leg given a settlement date

.. function:: has_occurred(cf::CashFlow, ref_date::Date, include_settlement_cf::Bool = true)

    Determines whether or not a particular cash flow has occurred or not

.. function:: maturity_date(leg::Leg)

    Returns the maturity date of a cash flow leg

.. function:: yield(leg::Leg, npv::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date, accuracy::Float64, max_iter::Int, guess::Float64)

    Calculates the yield of a cash flow leg


Zero Coupon Leg
~~~~~~~~~~~~~~~

.. code-block:: julia

    type ZeroCouponLeg <: Leg
      redemption::SimpleCashFlow
    end

.. function:: npv(leg::ZeroCouponLeg, yts::YieldTermStructure, settlement_date::Date, npv_date::Date)

    Calculates the npv of a zero coupon cash flow leg given a term structure, settlement date, and npv date

.. function:: npv(leg::ZeroCouponLeg, y::InterestRate, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date)

    Calculates the npv of a zero coupon cash flow leg given an interest rate, settlement date, and npv date

.. function:: duration(::ModifiedDuration, leg::ZeroCouponLeg, y::InterestRate, dc::DC, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date = Date())

    Calculates the modified duration of a zero coupon cash flow leg
