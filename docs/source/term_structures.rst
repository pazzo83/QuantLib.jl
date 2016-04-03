Term Structures and Curves
==========================

QuantLib.jl has various term structures and curves for asset pricing.

Various methods of bootstrapping rate curves are also available.


General Term Structure methods
------------------------------

.. function:: time_from_reference(ts::TermStructure, date::Date)

    Returns the time from the term structure's reference date given the day counter used by the term structure and a passed in date


Bootstrapping
-------------

QuantLib.jl has an iterative bootstrap type for bootstrapping a rate curve.

This bootstrapper uses a Brent Solver and Finite Differences Newton-Safe solver for bootstrap calculations

.. code-block:: julia

    type IterativeBootstrap <: Bootstrap
      firstSolver::BrentSolver
      solver::FiniteDifferenceNewtonSafe
    end

.. function:: IterativeBootstrap()

    Default constructor for the iterative bootstrap

.. function:: initialize(::IterativeBootstrap, ts::TermStructure)

    Initializes a term structure curve to prepare it for bootstrapping


Various traits govern how the bootstrapping is done, based on the curve that is being created.  The currently implemented traits are:

.. code-block:: julia

    type Discount <: BootstrapTrait end
    type HazardRate <: BootstrapTrait end


Bootstrap Helpers
~~~~~~~~~~~~~~~~~

QuantLib.jl's bootstrapping uses various helpers to construct the appropriate curve.


**FixedRateBondHelper**

Fixed-coupon bond helper for curve bootstrap

.. code-block:: julia

    type FixedRateBondHelper <: BondHelper
      price::Quote
      bond::FixedRateBond
    end


**SwapRateHelper**

Rate helper for bootstrapping over swap rates

.. code-block:: julia

    type SwapRateHelper <: RateHelper
      rate::Quote
      tenor::Dates.Period
      fwdStart::Dates.Period
      swap::VanillaSwap
    end


**DepositRateHelper**

Rate helper for bootstrapping over deposit rates

.. code-block:: julia

    type DepositRateHelper <: RateHelper
      rate::Quote
      tenor::TenorPeriod
      fixingDays::Int
      calendar::BusinessCalendar
      convention::BusinessDayConvention
      endOfMonth::Bool
      dc::DayCount
      iborIndex::IborIndex
      evaluationDate::Date
      referenceDate::Date
      earliestDate::Date
      maturityDate::Date
      fixingDate::Date
    end


**FraRateHelper**

Rate helper for bootstrapping over FRA rates

.. code-block:: julia

    type FraRateHelper <: RateHelper
      rate::Quote
      evaluationDate::Date
      periodToStart::Dates.Period
      iborIndex::IborIndex
      fixingDate::Date
      earliestDate::Date
      latestDate::Date
    end


**SpreadCDSHelper**

Spread-quoted CDS hazard rate bootstrap helper

.. code-block:: julia

    type SpreadCDSHelper <: AbstractCDSHelper
      runningSpread::Quote
      tenor::Dates.Period
      settlementDays::Int
      calendar::BusinessCalendar
      frequency::Frequency
      paymentConvention::BusinessDayConvention
      dc::DayCount
      recoveryRate::Float64
      schedule::Schedule
      discountCurve::YieldTermStructure
      settlesAccrual::Bool
      paysAtDefaultTime::Bool
      protectionStart::Date
      probability::AbstractDefaultProbabilityTermStructure
      swap::CreditDefaultSwap
    end

.. function:: SpreadCDSHelper(runningSpread::Quote, tenor::Dates.Period, settlementDays::Int, calendar::BusinessCalendar, frequency::Frequency, paymentConvention::BusinessDayConvention, rule::DateGenerationRule, dc::DayCount, recoveryRate::Float64, discountCurve::YieldTermStructure = NullYieldTermStructure(), settlesAccrual::Bool = true, paysAtDefaultTime::Bool = true)

    Constructor for the SpreadCDSHelper



Yield Term Structures
---------------------

Interest-rate term structure

.. function:: discount(yts::YieldTermStructure, date::Date)

    Returns the discount factor from a given date

.. function discount(yts::YieldTermStructure, time_frac::Float64)

    Returns the discount factor from a given period of time from the reference date

.. function:: zero_rate(yts::YieldTermStructure, date::Date, dc::DayCount, comp::CompoundingType, freq::Frequency = Annual())

    Returns the implied zero-yield rate for a given date (returns an InterestRate object)

.. function:: zero_rate(yts::YieldTermStructure, time_frac::Float64, comp::CompoundingType, freq::Frequency = Annual())

    Returns the implied zero-yield rate for a given time fraction from the reference date (returns an InterestRate object)

.. function:: forward_rate(yts::YieldTermStructure, date1::Date, date2::Date, dc::DayCount, comp::CompoundingType, freq::Frequency)

    Returns the forward interest rate between two dates

.. function:: forward_rate(yts::YieldTermStructure, date::Date, period::Integer, dc::DayCount, comp::CompoundingType, freq::Frequency)

    Returns the forward interest rate between a date and another date based on the passed-in time period

.. function:: forward_rate(yts::YieldTermStructure, time1::Float64, time2::Float64, comp::CompoundingType, freq::Frequency)

    Returns the forward interest rate between two time periods calculated from the reference date


FlatForwardTermStructure
~~~~~~~~~~~~~~~~~~~~~~~~

Flat interest-rate curve

.. code-block:: julia

    type FlatForwardTermStructure <: YieldTermStructure
      settlementDays::Int
      referenceDate::Date
      calendar::BusinessCalendar
      forward::Quote
      dc::DayCount
      comp::CompoundingType
      freq::Frequency
      rate::InterestRate
      jumpTimes::Vector{JumpTime}
      jumpDates::Vector{JumpDate}
    end

.. function:: FlatForwardTermStructure(settlement_days::Int, referenceDate::Date, calendar::BusinessCalendar, forward::Quote, dc::DayCount, comp::CompoundingType = ContinuousCompounding(), freq::Frequency = QuantLib.Time.Annual())

    Constructor for a FlatForwardTermStructure, with a quote used to generate the InterestRate object

.. function:: FlatForwardTermStructure(referenceDate::Date, calendar::BusinessCalendar, forward::Quote, dc::DayCount, comp::CompoundingType = ContinuousCompounding(), freq::Frequency = QuantLib.Time.Annual())

    Constructor for a FlatForwardTermStructure with no settlement days passed (defaults to 0) and a quote to generate the interest rate object

.. function:: FlatForwardTermStructure(settlementDays::Int, calendar::BusinessCalendar, forward::Quote, dc::DayCount, comp::CompoundingType = ContinuousCompounding(), freq::Frequency = QuantLib.Time.Annual())

    Constructor for a FlatForwardTermStructure with no reference date passed (will be calculated the first time it is requested) and a quote to generate the interest rate object

.. function:: FlatForwardTermStructure(referenceDate::Date, forward::Float64, dc::DayCount)

    Constructor for a FlatForwardTermStructure with only a reference date, forward rate, and day count passed in.  Will default with a TargetCalendar, ContinuousCompounding, and an Annual frequence

.. function:: FlatForwardTermStructure(referenceDate::Date, forward::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency)

    Constructor for a FlatForwardTermStructure with no settlement days or calendar passed (defaults to 0 and TargetCalendar, respectively)


InterpolatedDiscountCurve
~~~~~~~~~~~~~~~~~~~~~~~~~

YieldTermStructure based on interpolation of discount factors

.. code-block:: julia

    type InterpolatedDiscountCurve <: InterpolatedCurve
      settlementDays::Int
      referenceDate::Date
      dc::DayCount
      interp::Interpolation
      cal::BusinessCalendar
      dates::Vector{Date}
      times::Vector{Float64}
      data::Vector{Float64}
    end

.. function:: InterpolatedDiscountCurve(dates::Vector{Date}, discounts::Vector{Float64}, dc::DayCount, interpolator::Interpolation)

    Constructor for the InterpolatedDiscountCurve, with a given interpolation method


PiecewiseYieldCurve
~~~~~~~~~~~~~~~~~~~

Piecewise yield term structure.  This term structure is bootstrapped on a number of interest rate instruments which are passed as a vector of RateHelper instances. Their maturities mark the boundaries of the interpolated segments.

Each segment is determined sequentially starting from the earliest period to the latest and is chosen so that the instrument whose maturity marks the end of such segment is correctly repriced on the curve.

.. code-block:: julia

    type PiecewiseYieldCurve <: InterpolatedCurve{P, T}
      lazyMixin::LazyMixin
      settlementDays::Int
      referenceDate::Date
      instruments::Vector{BootstrapHelper}
      dc::DayCount
      interp::Interpolation
      trait::BootstrapTrait
      accuracy::Float64
      boot::Bootstrap
      times::Vector{Float64}
      dates::Vector{Date}
      data::Vector{Float64}
      errors::Vector{Function}
      validCurve::Bool
    end


FittedBondDiscountCurve
~~~~~~~~~~~~~~~~~~~~~~~

Discount curve fitted to a set of fixed-coupon bonds.

This class fits a discount function d(t) over a set of bonds, using a user defined fitting method. The discount function is fit in such a way so that all cashflows of all input bonds, when discounted using d(t), will reproduce the set of input bond prices in an optimized sense. Minimized price errors are weighted by the inverse of their respective bond duration.

The FittedBondDiscountCurve class acts as a generic wrapper, while its inner class FittingMethod provides the implementation details. Developers thus need only derive new fitting methods from the latter.


.. code-block:: julia

    type FittedBondDiscountCurve <: Curve
      lazyMixin::LazyMixin
      settlementDays::Int
      referenceDate::Date
      calendar::BusinessCalendar
      bonds::Vector{BondHelper}
      dc::DayCount
      fittingMethod::FittingMethod
      accuracy::Float64
      maxEvaluations::Int
      simplexLambda::Float64
    end


.. function:: FittedBondDiscountCurve(settlementDays::Int, referenceDate::Date, calendar::BusinessCalendar, bonds::Vector{BondHelper}, dc::DayCount, fittingMethod::FittingMethod, accuracy::Float64, maxEvaluations::Int, simplexLambda::Float64)

    Constructor for the FittedBondDiscountCurve


Fitting Methods
~~~~~~~~~~~~~~~

Common interface for all fitting methods:

.. code-block:: julia

    type FittingMethodCommons{T <: Real}
      solution::Vector{T}
      guessSolution::Vector{T}
      numberOfIterations::Int
      minimumCostValue::Float64
      weights::Vector{T}
      costFunction::FittingCost
    end


**ExponentialSplinesFitting**

Exponential-splines fitting method

.. code-block:: julia

    type ExponentialSplinesFitting <: FittingMethod
      constrainAtZero::Bool
      size::Int
      commons::FittingMethodCommons
    end

.. function:: ExponentialSplinesFitting(constrainAtZero::Bool, size::Int)

    Constructor for the ExponentialSplines fitting method


**SimplePolynomialFitting**

Simple polynomial fitting method

.. code-block:: julia

    type SimplePolynomialFitting <: FittingMethod
      constrainAtZero::Bool
      degree::Int
      size::Int
      commons::FittingMethodCommons
    end

.. function:: SimplePolynomialFitting(constrainAtZero::Bool, degree::Int, size::Int)

    Constructor for the SimplePolynomial fitting method


**NelsonSiegelFitting**

Nelson-Siegel fitting method

.. code-block:: julia

    type NelsonSiegelFitting <: FittingMethod
      constrainAtZero::Bool
      size::Int
      commons::FittingMethodCommons
    end

.. function:: NelsonSiegelFitting(size::Int)

    Constructor for the Nelson Siegel fitting method


**SvenssonFitting**

Svensson Fitting method

.. code-block:: julia

    type SvenssonFitting <: FittingMethod
      constrainAtZero::Bool
      size::Int
      commons::FittingMethodCommons
    end

.. function:: SvenssonFitting(size::Int)

    Constructor for the Svensson fitting method


**CubicBSplinesFitting**

CubicSpline B-splines fitting method

.. code-block:: julia

    type CubicBSplinesFitting <: FittingMethod
      constrainAtZero::Bool
      size::Int
      knots::Vector{Float64}
      splines::BSpline
      N::Int
      commons::FittingMethodCommons
    end

.. function:: CubicBSplinesFitting(constrainAtZero::Bool, knots::Vector{Float64}, size::Int)

    Default constructor for the Cubic BSplines fitting method



Credit Term Structures
----------------------

Term structures for calculating default probability in credit products


InterpolatedHazardRateCurve
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Default probability term structure based on interpolation of hazard rates

.. code-block:: julia

    type InterpolatedHazardRateCurve <: InterpolatedDefaultProbabilityCurve{P}
      settlementDays::Int
      referenceDate::Date
      dc::DayCount
      interp::Interpolation
      cal::BusinessCalendar
      dates::Vector{Date}
      times::Vector{Float64}
      data::Vector{Float64}
    end

.. function:: InterpolatedHazardRateCurve(dates::Vector{Date}, hazardRates::Vector{Float64}, dc::DayCount, interpolator::Interpolation)

    Constructor for the InterpolatedHazardRateCurve, given a set of dates and hazard rates


PiecewiseDefaultCurve
~~~~~~~~~~~~~~~~~~~~~

Piecewise default-probability term structure

This term structure is bootstrapped on a number of credit instruments which are passed as a vector of DefaultProbabilityHelper instances. Their maturities mark the boundaries of the interpolated segments.

Each segment is determined sequentially starting from the earliest period to the latest and is chosen so that the instrument whose maturity marks the end of such segment is correctly repriced on the curve.

.. code-block:: julia

    type PiecewiseDefaultCurve <: InterpolatedDefaultProbabilityCurve{P}
      lazyMixin::LazyMixin
      settlementDays::Int
      referenceDate::Date
      instruments::Vector{BootstrapHelper}
      dc::DayCount
      interp::Interpolation
      trait::BootstrapTrait
      accuracy::Float64
      boot::Bootstrap
      times::Vector{Float64}
      dates::Vector{Date}
      data::Vector{Float64}
      errors::Vector{Function}
      validCurve::Bool
    end

.. function:: PiecewiseDefaultCurve(referenceDate::Date, instruments::Vector{BootstrapHelper}, dc::DayCount, interp::Interpolation, trait::BootstrapTrait, accuracy::Float64, boot::Bootstrap = IterativeBootstrap())

    Constructor for a PiecewiseDefaultCurve



Volatility Term Structures
--------------------------

ConstantOptionVolatility
~~~~~~~~~~~~~~~~~~~~~~~~

Constant caplet volatility, no time-strike dependence

.. code-block:: julia

    type ConstantOptionVolatility <: OptionletVolatilityStructure
      settlementDays::Int
      referenceDate::Date
      calendar::BusinessCalendar
      bdc::BusinessDayConvention
      volatility::Float64
      dc::DayCount
    end

.. function:: ConstantOptionVolatility(settlementDays::Int, calendar::BusinessCalendar, bdc::BusinessDayConvention, volatility::Float64, dc::DayCount)

    Constructor for ConstantOptionVolatility, with floating reference date


ConstantSwaptionVolatility
~~~~~~~~~~~~~~~~~~~~~~~~~~

Constant swaption volatility, no time-strike dependence

.. code-block:: julia

    type ConstantSwaptionVolatility <: SwaptionVolatilityStructure
      settlementDays::Int
      referenceDate::Date
      calendar::BusinessCalendar
      bdc::BusinessDayConvention
      volatility::Quote
      dc::DayCount
    end

.. function:: ConstantSwaptionVolatility(settlementDays::Int, cal::BusinessCalendar, bdc::BusinessDayConvention, volatility::Quote, dc::DayCount)

    Constructor for ConstantSwaptionVolatility, with floating reference date


BlackConstantVol
~~~~~~~~~~~~~~~~

Constant Black volatility, no time-strike dependence

.. code-block:: julia

    type BlackConstantVol <: BlackVolTermStructure
      referenceDate::Date
      settlementDays::Int
      calendar::BusinessCalendar
      volatility::Quote
      dc::DayCount
    end

.. function:: BlackConstantVol(refDate::Date, cal::BusinessCalendar, volatility::Float64, dc::DayCount)

    Constructor for BlackConstantVol term structure, fixed reference date
