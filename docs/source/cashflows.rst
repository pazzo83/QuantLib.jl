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


Fixed Rate Coupon
~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type FixedRateCoupon{DC <: DayCount} <: Coupon
      couponMixin::CouponMixin{DC}
      paymentDate::Date
      nominal::Float64
      rate::InterestRate
    end

This is a coupon for fixed rate cash flows

.. function:: FixedRateCoupon{DC <: DayCount}(paymentDate::Date, faceAmount::Float64, rate::InterestRate, accrual_start::Date, accrual_end::Date, ref_start::Date, ref_end::Date, dc::DC, accrual::Float64)

    Constructor for FixedRateCoupon

.. function:: amount(coup::FixedRateCoupon)

    Calculates value of coupon

.. function:: calc_rate(coup::FixedRateCoupon) = coup.rate.rate

    Returns the coupon rate

.. function:: get_pay_dates(coups::Vector{Union{FixedRateCoupon, SimpleCashFlow}})

    Returns the pay dates as a vector of a vector of fixed rate coupons and simple cash flows (usually a redemption)

.. function:: get_reset_dates(coups::Vector{Union{FixedRateCoupon, SimpleCashFlow}})

    Returns the reset dates as a vector of a vector of fixed rate coupons and simple cash flows (usually a redemption)

.. function:: accrued_amount(coup::FixedRateCoupon, settlement_date::Date)

    Calculates the accrued amount of a fixed rate coupon given a settlement date


Ibor Coupon
~~~~~~~~~~~

.. code-block:: julia

    type IborCoupon{DC <: DayCount, X <: InterestRateIndex, ICP <: IborCouponPricer} <: Coupon
      couponMixin::CouponMixin{DC}
      paymentDate::Date
      nominal::Float64
      fixingDate::Date
      fixingValueDate::Date
      fixingEndDate::Date
      fixingDays::Int
      iborIndex::X
      gearing::Float64
      spread::Float64
      isInArrears::Bool
      spanningTime::Float64
      pricer::ICP
    end

.. function:: IborCoupon(paymentDate::Date, nominal::Float64, startDate::Date, endDate::Date, fixingDays::I, iborIndex::InterestRateIndex, gearing::Float64, spread::Float64, refPeriodStart::Date, refPeriodEnd::Date, dc::DayCount, isInArrears::Bool, pricer::IborCouponPricer)

    Constructor for IborCoupon

.. function:: amount(coup::IborCoupon)

    Calculates the Ibor coupon amount

.. function:: accrued_amount(coup::IborCoupon, settlement_date::Date)

    Calculates the accrued amount of the coupon given a settlement date

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


Fixed Rate Leg
~~~~~~~~~~~~~~

.. code-block:: julia

    type FixedRateLeg <: Leg
      coupons::Vector{Union{FixedRateCoupon, SimpleCashFlow}}
    end

.. function:: FixedRateLeg(schedule::Schedule, faceAmount::Float64, rate::Float64, calendar::BusinessCalendar, paymentConvention::BusinessDayConvention, dc::DayCount; add_redemption::Bool = true)

    Constructor for FixedRateLeg, passing in one rate

.. function:: FixedRateLeg(schedule::Schedule, faceAmount::Float64, rates::Vector{Float64}, calendar::BusinessCalendar, paymentConvention::BusinessDayConvention, dc::DayCount; add_redemption::Bool = false)

    Constructor for FixedRateleg, passing in a vector of rates


Ibor Leg
~~~~~~~~

.. code-block:: julia

    type IborLeg <: Leg
      coupons::Vector{Union{IborCoupon, SimpleCashFlow}}
    end

.. function:: IborLeg(schedule::Schedule, nominal::Float64, iborIndex::InterestRateIndex, paymentDC::DayCount, paymentAdj::BusinessDayConvention, fixingDays::Vector{Int} = fill(iborIndex.fixingDays, length(schedule.dates) - 1), gearings::Vector{Float64} = ones(length(schedule.dates) - 1), spreads::Vector{Float64} = zeros(length(schedule.dates) - 1), caps::Vector{Float64} = Vector{Float64}(), floors::Vector{Float64} = Vector{Float64}(), isInArrears::Bool = false, isZero::Bool = false, pricer::IborCouponPricer = BlackIborCouponPricer(); add_redemption::Bool = true)

    Constructor for Ibor Leg (will construct the coupons given the parameters passed in)
