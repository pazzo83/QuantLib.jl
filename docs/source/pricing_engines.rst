Pricing Engines
===============

Pricing engines are the main pricing tools in QuantLib.jl.  Each asset type has a variety of different pricing engines, depending on the pricing method.  Every asset is associated with a pricing engine , which is used to calculate NPV and other asset data.

Pricing engines usually have one or more term structures tied to them for pricing.

Bond Pricing Engines
--------------------

These pricing engines use a variety of methods to price bonds in QuantLib.jl

BlackCallableFixedRateBondEngine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Black-formula callable fixed rate bond engine

Callable fixed rate bond Black engine. The embedded (European) option follows the Black "European bond option" treatment in Hull, Fourth Edition, Chapter 20.

.. code-block:: julia

    type BlackCallableFixedRateBondEngine <: PricingEngine
      volatility::CallableBondVolatilityStructure
    end

.. function:: BlackCallableFixedRateBondEngine(fwdYieldVol::Quote)

    Constructor for the Black-Callable Fixed Rate Bond Engine.  This will construct the volatility term structure.


DiscountingBondEngine
~~~~~~~~~~~~~~~~~~~~~

Basic discounting bond engine using a yield term structure for pricing

.. code-block:: julia

    type DiscountingBondEngine <: PricingEngine
      yts::YieldTermStructure
    end

.. function:: DiscountingBondEngine()

    Optional constructor that will generate a dummy null term structure that must be changed later.


TreeCallableFixedRateEngine
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Numerical lattice engine for callable fixed rate bonds

.. code-block:: julia

    type TreeCallableFixedRateEngine <: LatticeShortRateModelEngine{S}
      model::ShortRateModel
      timeSteps::Int
      common::LatticeShortRateModelEngineCommon
    end

.. function:: TreeCallableFixedRateEngine(model::ShortRateModel, timeSteps::Int)

    Default constructor for the TreeCallableFixedRateEngine.  This will generate the necessary lattice for pricing.


Credit Pricing Engines
----------------------

These pricing engines are used for credit-related products

MidPointCdsEngine
~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type MidPointCdsEngine <: PricingEngine
      probability::AbstractDefaultProbabilityTermStructure
      recoveryRate::Float64
      discountCurve::YieldTermStructure
    end


Swap Pricing Engines
--------------------

Pricing engines used to price swaps


DiscountingSwapEngine
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type DiscountingSwapEngine <: PricingEngine
      yts::YieldTermStructure
      includeSettlementDateFlows::Bool
    end

.. function:: DiscountingSwapEngine(yts::YieldTermStructure, includeSettlementDateFlows::Bool = true)

    Constructor for the DiscountingSwapEngine

.. function:: DiscountingSwapEngine()

    Constructor for the DiscountingSwapEngine that will generate a dummy null yield term structure.  This must be changed for any pricing.


Swaption Pricing Engines
------------------------

Pricing engines to price swaptions


BlackSwaptionEngine
~~~~~~~~~~~~~~~~~~~

Shifted Lognormal Black-formula swaption engine.  The engine assumes that the exercise date equals the start date of the passed swap.

.. code-block:: julia

    type BlackSwaptionEngine <: PricingEngine
      yts::YieldTermStructure
      vol::Quote
      volStructure::SwaptionVolatilityStructure
      dc::DayCount
      displacement::Float64
    end

.. function:: BlackSwaptionEngine(yts::YieldTermStructure, vol::Quote, dc::DayCount, displacement::Float64 = 0.0)

    Constructor for the BlackSwaptionEngine.  Will generate the volatility term structure.


FdHullWhiteSwaptionEngine
~~~~~~~~~~~~~~~~~~~~~~~~~

Finite Differences swaption engine using a Hull White Model

.. code-block:: julia

    type FdHullWhiteSwaptionEngine <: PricingEngine
      model::HullWhite
      tGrid::Int
      xGrid::Int
      dampingSteps::Int
      invEps::Float64
      schemeDesc::FdmSchemeDesc
      ts::TermStructure
    end

.. function:: FdHullWhiteSwaptionEngine(model::HullWhite, tGrid::Int = 100, xGrid::Int = 100, dampingSteps::Int = 0, invEps::Float64 = 1e-5, schemeDesc::FdmSchemeDesc = FdmSchemeDesc(Douglas()))

    Constructor for the FdHullWhiteSwaptionEngine.  Uses the term structure from the hull white model by default.


FdG2SwaptionEngine
~~~~~~~~~~~~~~~~~~

Finite Differences swaption engine using a two-factor Gaussian model (G2)

.. code-block:: julia

    type FdG2SwaptionEngine <: PricingEngine
      model::G2
      tGrid::Int
      xGrid::Int
      yGrid::Int
      dampingSteps::Int
      invEps::Float64
      schemeDesc::FdmSchemeDesc
      ts::TermStructure
    end

.. function:: FdG2SwaptionEngine(model::G2, tGrid::Int = 100, xGrid::Int = 50, yGrid::Int = 50, dampingSteps::Int = 0, invEps::Float64 = 1e-5, schemeDesc::FdmSchemeDesc = FdmSchemeDesc(Hundsdorfer()))

    Constructor for the FdG2SwaptionEngine.  Uses the term structure from the G2 model by default.


G2SwaptionEngine
~~~~~~~~~~~~~~~~

Swaption priced by means of the Black formula, using a G2 model.  The engine assumes that the exercise date equals the start date of the passed swap.

.. code-block:: julia

    type G2SwaptionEngine <: PricingEngine
      model::G2
      range::Float64
      intervals::Int
    end


Gaussian1DSwaptionEngine
~~~~~~~~~~~~~~~~~~~~~~~~

One factor gaussian model swaption engine.  All fixed coupons with start date greater or equal to the respective option expiry are considered to be part of the exercise into right.

Cash settled swaptions are not supported.

.. code-block:: julia

    abstract GaussianProbabilities
    type NoneProbabilities <: GaussianProbabilities end
    type NaiveProbabilities <: GaussianProbabilities end
    type DigitalProbabilities <: GaussianProbabilities end

    type Gaussian1DSwaptionEngine <: PricingEngine
      model::Gaussian1DModel
      integrationPoints::Int
      stddevs::Float64
      extrapolatePayoff::Bool
      flatPayoffExtrapolation::Bool
      discountCurve::YieldTermStructure
      probabilities::GaussianProbabilities
    end

.. function:: Gaussian1DSwaptionEngine(model::Gaussian1DModel, integrationPoints::Int = 64, stddevs::Float64 = 7.0, extrapolatePayoff::Bool = true, flatPayoffExtrapolation::Bool = false, discountCurve::YieldTermStructure = NullYieldTermStructure(), probabilities::GaussianProbabilities = NoneProbabilities())

    Constructor for the Gaussian 1D Swaption Engine.


Gaussian1DNonstandardSwaptionEngine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

One factor gaussian model non standard swaption engine.

All fixed coupons with start date greater or equal to the respective option expiry are considered to be part of the exercise into right.

All float coupons with start date greater or equal to the respective option expiry are considered to be part of the exercise into right.

For redemption flows an associated start date is considered in the criterion, which is the start date of the regular xcoupon period with same payment date as the redemption flow.

Cash settled swaptions are not supported

.. code-block:: julia

    type Gaussian1DNonstandardSwaptionEngine <: PricingEngine
      model::Gaussian1DModel
      integrationPoints::Int
      stddevs::Float64
      extrapolatePayoff::Bool
      flatPayoffExtrapolation::Bool
      oas::Quote
      discountCurve::YieldTermStructure
      probabilities::GaussianProbabilities
    end

.. function:: Gaussian1DNonstandardSwaptionEngine(model::Gaussian1DModel, integrationPoints::Int = 64, stddevs::Float64 = 7.0, extrapolatePayoff::Bool = true, flatPayoffExtrapolation::Bool = false, oas::Quote = Quote(-1.0), discountCurve::YieldTermStructure = NullYieldTermStructure(), probabilities::GaussianProbabilities = NoneProbabilities())

    Constructor for the Gaussian 1D Non-standard Swaption Engine.


JamshidianSwaptionEngine
~~~~~~~~~~~~~~~~~~~~~~~~

Swaption engine using Jamshidian's decomposition.  Concerning the start delay cf. `<http://ssrn.com/abstract=2246054>`_

.. code-block:: julia

    type JamshidianSwaptionEngine <: PricingEngine
      model::ShortRateModel
    end


TreeSwaptionEngine
~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type TreeSwaptionEngine <: LatticeShortRateModelEngine{S}
      model::ShortRateModel
      timeSteps::Int
      common::LatticeShortRateModelEngineCommon
    end

.. function:: TreeSwaptionEngine(model::ShortRateModel, tg::TimeGrid)

    Constructor for the TreeSwaptionEngine, using a time grid.  This will generate the necessary lattice from the time grid.

.. function:: TreeSwaptionEngine(model::ShortRateModel, timeSteps::Int)

    Constructor for the TreeSwaptionEngine, using a number of time steps.  With this construction, the necessary tree will not be generated until calculation.


Vanilla Option Pricing Engines
------------------------------

QuantLib.jl has a number of methods to price European, American, and Bermudan options


AnalyticEuropeanEngine
~~~~~~~~~~~~~~~~~~~~~~

Pricing engine for European vanilla options using analytical formulae

.. code-block:: julia

    type AnalyticEuropeanEngine <: PricingEngine
      process::AbstractBlackScholesProcess
    end


AnalyticHestonEngine
~~~~~~~~~~~~~~~~~~~~

Heston-model engine based on Fourier transform

.. code-block:: julia

    type AnalyticHestonEngine <: AbstractHestonEngine
      model::HestonModel
      evaluations::Int
      cpxLog::ComplexLogFormula
      integration::HestonIntegration
    end

.. function:: AnalyticHestonEngine(hestonModel::HestonModel)

    Constructor for the AnalyticHestonEngine.  Will construct the integration method as well.


BaroneAdesiWhaleyApproximationEngine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Barone-Adesi and Whaley pricing engine for American options (1987)

.. code-block:: julia

    type BaroneAdesiWhaleyApproximationEngine <: PricingEngine
      process::AbstractBlackScholesProcess
    end


BatesEngine
~~~~~~~~~~~

Bates model engines based on Fourier transform

.. code-block:: julia

    type BatesEngine <: AbstractHestonEngine
      model::BatesModel
      evaluations::Int
      cpxLog::ComplexLogFormula
      integration::HestonIntegration
    end

.. function:: BatesEngine(batesModel::BatesModel)

    Constructor for the BatesEngine.  This will also construct the integration method.


BinomialVanillaEngine
~~~~~~~~~~~~~~~~~~~~~

Pricing engine for vanilla options using binomial trees

.. code-block:: julia

    type BinomialVanillaEngine{P <: AbstractBlackScholesProcess, T <: BinomialTreeType} <: AbstractVanillaEngine
      process::P
      timeSteps::Int
      treeClass::Type{T}
    end


BjerksundStenslandApproximationEngine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Bjerksund and Stensland pricing engine for American options (1993)

.. code-block:: julia

    type BjerksundStenslandApproximationEngine <: PricingEngine
      process::AbstractBlackScholesProcess
    end


FDEuropeanEngine
~~~~~~~~~~~~~~~~

Pricing engine for European options using finite-differences

.. code-block:: julia

    type FDEuropeanEngine <: AbstractFDVanillaEngine
      process::AbstractBlackScholesProcess
      timeSteps::Int
      gridPoints::Int
      timeDependent::Bool
      exerciseDate::Date
      finiteDifferenceOperator::TridiagonalOperator
      intrinsicValues::SampledCurve
      BCs::Vector{BoundaryCondition}
      sMin::Float64
      center::Float64
      sMax::Float64
      prices::SampledCurve
      fdEvolverFunc::Function
    end

.. function:: FDEuropeanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100, gridPoints::Int = 100, timeDependent::Bool = false)

    Constructor for the FDEuropeanEngine, requires a finite-differences evolver


FDBermudanEngine
~~~~~~~~~~~~~~~~

Finite-differences Bermudan engine

.. code-block:: julia

    type FDBermudanEngine <: FDMultiPeriodEngine
      process::AbstractBlackScholesProcess
      timeSteps::Int
      gridPoints::Int
      timeDependent::Bool
      exerciseDate::Date
      finiteDifferenceOperator::TridiagonalOperator
      intrinsicValues::SampledCurve
      BCs::Vector{BoundaryCondition}
      sMin::Float64
      center::Float64
      sMax::Float64
      prices::SampledCurve
      stoppingTimes::Vector{Float64}
      timeStepPerPeriod::Int
      fdEvolverFunc::Function
    end

.. function:: FDBermudanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100, gridPoints::Int = 100, timeDependent::Bool = false)

    Constructor for the FDBermudanEngine, requires a finite-differences evolver


FDAmericanEngine
~~~~~~~~~~~~~~~~

Finite-differences pricing engine for American one asset options

.. code-block:: julia

    type FDAmericanEngine <: FDStepConditionEngine
      process::AbstractBlackScholesProcess
      timeSteps::Int
      gridPoints::Int
      timeDependent::Bool
      exerciseDate::Date
      finiteDifferenceOperator::TridiagonalOperator
      intrinsicValues::SampledCurve
      BCs::Vector{BoundaryCondition}
      sMin::Float64
      center::Float64
      sMax::Float64
      prices::SampledCurve
      controlPrices::SampledCurve
      controlBCs::Vector{BoundaryCondition}
      controlOperator::TridiagonalOperator
      fdEvolverFunc::Function
    end


.. function:: FDAmericanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100, gridPoints::Int = 100, timeDependent::Bool = false)

    Constructor for the FDAmericanEngine, requires a finite-differences evolver


IntegralEngine
~~~~~~~~~~~~~~

Pricing engine for European vanilla options using integral approach

.. code-block:: julia

    type IntegralEngine <: PricingEngine
      process::AbstractBlackScholesProcess
    end


MCAmericanEngine
~~~~~~~~~~~~~~~~

American Monte Carlo engine, using the Longstaff-Schwarz Monte Carlo engine for early exercise options

.. code-block:: julia

    type MCAmericanEngine <: MCLongstaffSchwartzEngine
      process::AbstractBlackScholesProcess
      timeSteps::Int
      timeStepsPerYear::Int
      requiredSamples::Int
      maxSamples::Int
      requiredTolerance::Float64
      brownianBridge::Bool
      seed::Int
      nCalibrationSamples::Int
      polynomOrder::Int
      polynomType::LsmBasisSystemPolynomType
      mcSimulation::MCSimulation
      pathPricer::LongstaffSchwartzPathPricer
    end

.. function:: MCAmericanEngine(process::AbstractBlackScholesProcess; timeSteps::Int = -1, timeStepsPerYear::Int = -1, brownianBridge::Bool = false, antitheticVariate::Bool = false, requiredSamples::Int = -1, requiredTolerance::Float64 = -1.0, maxSamples::Int = typemax(Int), seed::Int = 0, rsg::AbstractRandomSequenceGenerator = InverseCumulativeRSG(seed), nCalibrationSamples::Int = 2048, polynomOrder::Int = 2, polynomType::LsmBasisSystemPolynomType = Monomial())

    Constructor for the MCAmericanEngine


MCEuropeanEngine
~~~~~~~~~~~~~~~~

European option pricing engine using Monte Carlo simulation

.. code-block:: julia

    type MCEuropeanEngine <: MCVanillaEngine
      process::AbstractBlackScholesProcess
      timeSteps::Int
      timeStepsPerYear::Int
      requiredSamples::Int
      maxSamples::Int
      requiredTolerance::Float64
      brownianBridge::Bool
      seed::Int
      mcSimulation::MCSimulation
    end

.. function:: MCEuropeanEngine(process::AbstractBlackScholesProcess; timeSteps::Int = -1, timeStepsPerYear::Int = -1, brownianBridge::Bool = false, antitheticVariate::Bool = false, requiredSamples::Int = -1, requiredTolerance::Float64 = -1.0, maxSamples::Int = typemax(Int), seed::Int = 1, rsg::AbstractRandomSequenceGenerator = InverseCumulativeRSG(seed))

    Constructor for the MCEuropeanEngine


Black Calculator
----------------

Black 1976 calculator type

.. code-block:: julia

    type BlackCalculator
      payoff::StrikedTypePayoff
      strike::Float64
      forward::Float64
      stdDev::Float64
      discount::Float64
      variance::Float64
      d1::Float64
      d2::Float64
      alpha::Float64
      beta::Float64
      DalphaDd1::Float64
      DbetaDd2::Float64
      n_d1::Float64
      cum_d1::Float64
      n_d2::Float64
      cum_d2::Float64
      x::Float64
      DxDs::Float64
      DxDstrike::Float64
    end

.. function:: BlackCalculator(p::StrikedTypePayoff, fwd::Float64, stdDev::Float64, disc::Float64)

    Constructor for the Black Calculator

.. function:: initialize!(calc::BlackCalculator, p::StrikedTypePayoff)

    Initializes the black calculator with a given payoff

.. function:: value(calc::BlackCalculator)

    Returns value of black calculator

.. function:: delta(calc::BlackCalculator, spot::Float64)

    Sensitivity to change in the underlying spot price.

.. function:: vega(calc::BlackCalculator, mat::Float64)

    Sensitivity to volatility.


Discretized Assets
------------------

Discretized Assets used by numerical methods


DiscretizedSwap
~~~~~~~~~~~~~~~

.. code-block:: julia

    type DiscretizedSwap <: DiscretizedAsset
      nominal::Float64
      swapT::SwapType
      fixedResetTimes::Vector{Float64}
      fixedPayTimes::Vector{Float64}
      floatingResetTimes::Vector{Float64}
      floatingPayTimes::Vector{Float64}
      args::VanillaSwapArgs
      common::DiscretizedAssetCommon
    end

.. function:: DiscretizedSwap(nominal::Float64, swapT::SwapType, referenceDate::Date, dc::DayCount, fixedPayDates::Vector{Date}, fixedResetDates::Vector{Date}, floatingPayDates::Vector{Date}, floatingResetDates::Vector{Date}, args::VanillaSwapArgs)

    Constructor for a Discretized Swap


DiscretizedSwaption
~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type DiscretizedSwaption <: DiscretizedOption
      underlying::DiscretizedSwap
      exercise::Exercise
      exerciseTimes::Vector{Float64}
      fixedPayDates::Vector{Date}
      fixedResetDates::Vector{Date}
      floatingPayDates::Vector{Date}
      floatingResetDates::Vector{Date}
      lastPayment::Float64
      common::DiscretizedAssetCommon
    end

.. function:: DiscretizedSwaption(swaption::Swaption, referenceDate::Date, dc::DayCount)

    Constructor for a Discretized Swaption, based on an underlying swaption
