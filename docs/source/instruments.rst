Instruments
===========

These are some of the core types of QuantLib.jl.  They involve the various instruments that QuantLib.jl is used to price, such as bonds, swaps, and options.

.. code-block:: julia

    abstract Instrument <: LazyObject

Bonds
-----

.. code-block:: julia

    abstract Bond <: Instrument

Common Bond type and functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type BondMixin
        settlementDays::Int
        issueDate::Date
        maturityDate::Date
    end

The bond mixin provides a common interface for all bonds.  These attributes are shared by all the bond types in QuantLib.jl

.. function:: get_settlement_days(bond::Bond)

    Returns the bond settlement days

.. function:: get_issue_date(bond::Bond)

    Returns the issue date of the bond

.. function:: get_maturity_date(bond::Bond)

    Returns the bond's maturity date

.. function:: get_settlement_date(b::Bond)

    Returns the bond's settlement date

.. function:: notional(bond::Bond, d::Date)

    Returns the bond's notional

.. function:: accrued_amount(bond::Bond, settlement::Date)

    Calculates the accrued amount of a bond's cashflows given a settlement date

.. function:: maturity_date(bond::Bond)

    Returns the maturity date of the bond

.. function:: yield(bond::Bond, clean_price::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, settlement::Date, accuracy::Float64 = 1.0e-10, max_iter::Int = 100, guess::Float64 = 0.05)

    Calculates the bond yield, passing in an accuracy limit and iteration limit for the calculation

.. function:: yield(bond::Bond, dc::DayCount, compounding::CompoundingType, freq::Frequency, settlement::Date, accuracy::Float64 = 1.0e-10, max_iter::Int = 100)

    Calculates the bond yield as above, but with no clean price passed in (this gets calculated)

.. function:: yield(bond::Bond, dc::DayCount, compounding::CompoundingType, freq::Frequency, accuracy::Float64 = 1.0e-10, max_iter::Int = 100)

    Calculates the bond yield as above, but with no clean price or settlement date passed in

.. function:: duration(bond::Bond, yld::InterestRate, duration_::Duration, dc::DayCount, settlement_date::Date)

    Calculates the bond duration

.. function:: duration(bond::Bond, yld::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, duration_::Duration, settlement_date::Date)

    Calculates the bond duration, with a rate passed in instead of an InterestRate object

.. function:: npv(bond::Bond)

    Calculates the bond NPV (this will trigger the bond calculation)

.. function:: clean_price(bond::Bond, settlement_value::Float64, settlement_date::Date)

    Calculates the bond clean price given a settlement value and settlement date

.. function:: clean_price(bond::Bond)

    Calculates the bond clean price (this will trigger calculation to get the settlement value)

.. function:: dirty_price(bond::Bond, settlement_value::Float64, settlement_date::Date)

    Calculates the bond dirty price given a settlement value and settlement date

.. function:: dirty_price(bond::Bond)

    Calculates the bond dirty price

.. function:: settlement_date(bond::Bond, d::Date = Date())

    Returns the bond settlement date


Fixed Rate Bond
~~~~~~~~~~~~~~~

.. code-block:: julia

    type FixedRateBond <: Bond
      lazyMixin::LazyMixin
      bondMixin::BondMixin
      faceAmount::Float64
      schedule::Schedule
      cashflows::FixedRateLeg
      dc::DayCount
      redemption::Float64
      startDate::Date
      pricingEngine::PricingEngine
      settlementValue::Float64
    end

.. function:: FixedRateBond(settlementDays::Int, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DayCount, paymentConvention::BusinessDayConvention, redemption::Float64, issueDate::Date, calendar::BusinessCalendar, pricing_engine::PricingEngine)

    Constructor for a Fixed Rate Bond


Floating Rate Bond
~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type FloatingRateBond <: Bond
      lazyMixin::LazyMixin
      bondMixin::BondMixin
      faceAmount::Float64
      schedule::Schedule
      cashflows::IborLeg
      iborIndex::InterestRateIndex
      dc::DayCount
      fixingDays::Int
      gearings::Vector{Float64}
      spreads::Vector{Float64}
      caps::Vector{Float64}
      floors::Vector{Float64}
      inArrears::Bool
      redemption::Float64
      pricingEngine::PricingEngine
      settlementValue::Float64
    end

.. function:: FloatingRateBond(settlementDays::Int, faceAmount::Float64, schedule::Schedule, iborIndex::InterestRateIndex, dc::DayCount, convention::BusinessDayConvention, fixingDays::Int, issueDate::Date, pricingEngine::PricingEngine, inArrears::Bool = false, redemption::Float64 = 100.0, gearings::Vector{Float64} = ones(length(schedule.dates) - 1), spreads::Vector{Float64} = zeros(length(schedule.dates) - 1), caps::Vector{Float64} = Vector{Float64}(), floors::Vector{Float64} = Vector{Float64}(); cap_vol::OptionletVolatilityStructure=NullOptionletVolatilityStructure())

    Constructor for a floating rate bond, optional argument to pass in an optionlet volatility term structure for the floating rate coupon pricer


Zero Coupon Bond
~~~~~~~~~~~~~~~~

.. code-block:: julia

    type ZeroCouponBond <: Bond
      lazyMixin::LazyMixin
      bondMixin::BondMixin
      faceAmount::Float64
      redemption::Float64
      cashflows::ZeroCouponLeg
      calendar::BusinessCalendar
      settlementValue::Float64
      pricingEngine::PricingEngine
    end


.. function:: ZeroCouponBond(settlementDays::Int, calendar::BusinessCalendar, faceAmount::Float64, maturityDate::Date, paymentConvention::BusinessDayConvention=Following(), redemption::Float64=100.0, issueDate::Date=Date(), pe::PricingEngine = DiscountingBondEngine())

    Constructor for a zero coupon bond


Callable Fixed Rate Bond
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    abstract AbstractCallableBond <: Bond

    type CallableFixedRateBond <: AbstractCallableBond
      lazyMixin::LazyMixin
      bondMixin::BondMixin
      faceAmount::Float64
      schedule::Schedule
      cashflows::FixedRateLeg
      dc::DayCount
      redemption::Float64
      startDate::Date
      pricingEngine::PricingEngine
      settlementValue::Float64
      putCallSchedule::CallabilitySchedule
      blackEngine::PricingEngine
      blackVolQuote::Quote
    end

.. function:: CallableFixedRateBond(settlementDays::Int, faceAmount::Float64, schedule::Schedule, coupons::Union{Vector{Float64}, Float64}, accrualDayCounter::DayCount, paymentConvention::BusinessDayConvention, redemption::Float64, issueDate::Date, putCallSchedule::CallabilitySchedule, pe::PricingEngine)

    Constructor for a callable fixed rate bond


Callability Schedule
~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type DirtyCall <: CallType end
    type CleanCall <: CallType end

    type Price
      amount::Float64
      callType::CallType
    end

    type Callability
      price::Price
      optionType::OptionType
      date::Date
    end

    typealias CallabilitySchedule Vector{Callability}



Forwards
--------

.. code-block:: julia

    abstract AbstractForward <: Instrument

    type ForwardRateAgreement <: AbstractForward
      lazyMixin::LazyMixin
      underlyingIncome::Float64
      underlyingSpotValue::Float64
      dc::DayCount
      calendar::BusinessCalendar
      convention::BusinessDayConvention
      settlementDays::Int
      payoff::ForwardTypePayoff
      valueDate::Date
      maturityDate::Date
      discountCurve::YieldTermStructure
      fraType::PositionType
      forwardRate::InterestRate
      strikeForwardRate::InterestRate
      notionalAmount::Float64
      iborIndex::IborIndex
    end


.. function:: ForwardRateAgreement(valueDate::Date, maturityDate::Date, position::PositionType, strikeForward::Float64, notionalAmount::Float64, iborIndex::IborIndex, discountCurve::YieldTermStructure)

    Constructor for a forward rate agreement

.. function:: spot_value(fra::ForwardRateAgreement)

    Calculates the spot value of a forward rate agreement (triggers calculation)

.. function:: forward_value(fra::ForwardRateAgreement)

    Calculates the forward value of a forward rate agreement (triggers calculation)

.. function:: forward_rate(fra::ForwardRateAgreement)

    Calculates the forward rate of a forward rate agreement (triggers calculation)

.. function:: implied_yield(fra::ForwardRateAgreement, underlyingSpotValue::Float64, forwardValue::Float64, settlementDate::Date, comp::CompoundingType, dc::DayCount)

    Calculates the implied yield of a forward rate agreement (triggers calculation)


Options
-------

.. code-block:: julia

    abstract Option{E} <: Instrument
    abstract OneAssetOption{E} <: Option{E}

    typealias EuropeanOption Option{EuropeanExercise}
    typealias AmericanOption Option{AmericanExercise}
    typealias BermudanOption Option{BermudanExercise}

    type VanillaOption <: OneAssetOption{Exercise}
      lazyMixin::LazyMixin
      payoff::StrikedTypePayoff
      exercise::Exercise
      pricingEngine::PricingEngine
      results::OptionResults
    end

    # These are the results that are calculated when the option is priced
    type OptionResults # Greeks
      delta::Float64
      gamma::Float64
      theta::Float64
      vega::Float64
      rho::Float64
      dividendRho::Float64
      deltaForward::Float64
      elasticity::Float64
      thetaPerDay::Float64
      strikeSensitivity::Float64
      itmCashProbability::Float64
      value::Float64
    end

.. function:: VanillaOption(payoff::StrikedTypePayoff, exercise::Exercise, pe::PricingEngine)

    Constructor for a vanilla option


Payoffs
-------

.. code-block:: julia

    abstract AbstractPayoff
    abstract StrikedTypePayoff <: AbstractPayoff

    type Put <: OptionType end
    type Call <: OptionType end

    type PlainVanillaPayoff <: StrikedTypePayoff
      optionType::OptionType
      strike::Float64
    end


Swaps
-----

.. code-block:: julia

    abstract Swap <: Instrument

Vanilla Swaps
~~~~~~~~~~~~~

.. code-block:: julia

    type Payer <: SwapType end
    type Receiver <: SwapType end

    # these are the results that are calculated when the swap is priced
    type SwapResults <: Results
      legNPV::Vector{Float64}
      legBPS::Vector{Float64}
      npvDateDiscount::Float64
      startDiscounts::Vector{Float64}
      endDiscounts::Vector{Float64}
      fairRate::Float64
      value::Float64
    end

    type VanillaSwap <: Swap
      lazyMixin::LazyMixin
      swapT::SwapType
      nominal::Float64
      fixedSchedule::Schedule
      fixedRate::Float64
      fixedDayCount::DayCount
      iborIndex::IborIndex
      spread::Float64
      floatSchedule::Schedule
      floatDayCount::DayCount
      paymentConvention::BusinessDayConvention
      legs::Vector{Leg}
      payer::Vector{Float64}
      pricingEngine::PricingEngine
      results::SwapResults
      args::VanillaSwapArgs
    end

.. function:: VanillaSwap(swapT::SwapType, nominal::Float64, fixedSchedule::Schedule, fixedRate::Float64, fixedDayCount::DayCount, iborIndex::IborIndex, spread::Float64, floatSchedule::Schedule, floatDayCount::DayCount, pricingEngine::PricingEngine = NullPricingEngine(), paymentConvention::B = floatSchedule.convention)

    Constructor for a vanilla swap

.. function:: fair_rate(swap::VanillaSwap)

    Calculates the fair rate of a vanilla swap (triggers calculation)


Nonstandard Swap
~~~~~~~~~~~~~~~~

This swap is used in all Gaussian Short Rate model calculations

.. code-block:: julia

    type NonstandardSwap <: Swap
        lazyMixin::LazyMixin
        swapT::SwapType
        fixedNominal::Vector{Float64}
        floatingNominal::Vector{Float64}
        fixedSchedule::Schedule
        fixedRate::Vector{Float64}
        fixedDayCount::DayCount
        iborIndex::IborIndex
        spread::Float64
        gearing::Float64
        floatSchedule::Schedule
        floatDayCount::DayCount
        paymentConvention::BusinessDayConvention
        intermediateCapitalExchange::Bool
        finalCapitalExchange::Bool
        legs::Vector{Leg}
        payer::Vector{Float64}
        pricingEngine::PricingEngine
        results::SwapResults
        args::NonstandardSwapArgs
    end

.. function:: NonstandardSwap(vs::VanillaSwap)

    Constructor for a Nonstandard Swap


Credit Default Swap
~~~~~~~~~~~~~~~~~~~

.. code-block:: julia
    type Buyer <: CDSProtectionSide end
    type Seller <: CDSProtectionSide end

    # These are calculated when the CDS is priced
    type CDSResults
      upfrontNPV::Float64
      couponLegNPV::Float64
      defaultLegNPV::Float64
      fairSpread::Float64
      fairUpfront::Float64
      couponLegBPS::Float64
      upfrontBPS::Float64
      value::Float64
    end

    type CreditDefaultSwap <: Swap
      lazyMixin::LazyMixin
      side::CDSProtectionSide
      notional::Float64
      spread::Float64
      schedule::Schedule
      convention::BusinessDayConvention
      dc::DayCount
      leg::FixedRateLeg
      upfrontPayment::SimpleCashFlow
      settlesAccrual::Bool
      paysAtDefaultTime::Bool
      protectionStart::Date
      pricingEngine::PricingEngine
      claim::FaceValueClaim
      results::CDSResults
    end

.. function:: CreditDefaultSwap(side::CDSProtectionSide, notional::Float64, spread::Float64, schedule::Schedule, convention::BusinessDayConvention, dc::DayCount, settlesAccrual::Bool, paysAtDefaultTime::Bool, protectionStart::Date, pricingEngine::PricingEngine)

    Constructor for a credit default swap


Swaptions
---------

.. code-block:: julia

    type SettlementPhysical <: SettlementType end
    type SettlementCash <: SettlementType end

    type SwaptionResults{S <: AbstractString}
      value::Float64
      additionalResults::Dict{S, Float64}
    end

    type Swaption <: Option
      lazyMixin::LazyMixin
      swap::VanillaSwap
      exercise::Exercise
      delivery::SettlementType
      results::SwaptionResults
      pricingEngine::PricingEngine
    end

.. function:: Swaption(swap::VanillaSwap, exercise::Exercise)

    Constructor for a swaption with no pricing engine or settlement type set

.. function:: Swaption(swap::VanillaSwap, exercise::Exercise, pe::PricingEngine)

    Constructor for a swaption with no settlement type set


Nonstandard Swaptions
~~~~~~~~~~~~~~~~~~~~~

These swaptions are used in Gaussian Short Rate model calculations

.. code-block:: julia

    type NonstandardSwaption <: Option
      lazyMixin::LazyMixin
      swap::NonstandardSwap
      exercise::Exercise
      delivery::SettlementType
      results::SwaptionResults
      pricingEngine::PricingEngine
    end

.. function:: NonstandardSwaption(swap::NonstandardSwap, exercise::Exercise)

    Constructor for a nonstandard swaption with no pricing engine or settlement type set

.. function:: NonstandardSwaption(swap::NonstandardSwap, exercise::Exercise, pe::PricingEngine)

    Constructor for a nonstandard swaption with no settlement type set

.. function calibration_basket(swaption::NonstandardSwaption, swaptionEngine::PricingEngine, swapIndex::SwapIndex, swaptionVol::SwaptionVolatilityStructure, basketType::CalibrationBasketType)

    Builds a basket for model calibration (calls the swaption engine's calibration_basket method)
