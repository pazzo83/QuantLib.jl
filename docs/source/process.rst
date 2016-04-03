Stochastic Processes
====================

Stochastic processes for use in various models and pricing engines


Heston Processes
----------------

HestonProcess
~~~~~~~~~~~~~

Square-root stochastic-volatility Heston process

.. code-block:: julia

    type HestonProcess <: AbstractHestonProcess
      s0::Quote
      riskFreeRate::YieldTermStructure
      dividendYield::YieldTermStructure
      v0::Float64
      kappa::Float64
      theta::Float64
      sigma::Float64
      rho::Float64
      disc::AbstractDiscretization
      hestonDisc::AbstractHestonDiscretization
    end

.. function:: HestonProcess(riskFreeRate::YieldTermStructure, dividendYield::YieldTermStructure, s0::Quote, v0::Float64, kappa::Float64, theta::Float64, sigma::Float64, rho::Float64, d::AbstractHestonDiscretization = QuadraticExponentialMartingale())

    Default constructor for the Heston Process


BatesProcess
~~~~~~~~~~~~

Square-root stochastic-volatility Bates process

.. code-block:: julia

    type BatesProcess<: AbstractHestonProcess
      hestonProcess::HestonProcess
      lambda::Float64
      nu::Float64
      delta::Float64
      m::Float64
    end

.. function:: BatesProcess(riskFreeRate::YieldTermStructure, dividendYield::YieldTermStructure, s0::Quote, v0::Float64, kappa::Float64, theta::Float64, sigma::Float64, rho::Float64, lambda::Float64, nu::Float64, delta::Float64, d::AbstractHestonDiscretization = FullTruncation())

    Constructor for the Bates Process.  Builds underlying Heston process as well


Black Scholes Processes
-----------------------

Various types of Black-Scholes stochastic processes

Black Scholes types are based off of a GeneralizedBlackScholesProcess, with a structure seen here:

.. code-block:: julia

    type GeneralizedBlackScholesProcess <: AbstractBlackScholesProcess
      x0::Quote
      riskFreeRate::YieldTermStructure
      dividendYield::YieldTermStructure
      blackVolatility::BlackVolTermStructure
      localVolatility::LocalConstantVol
      disc::AbstractDiscretization
      isStrikeDependent::Bool
      blackScholesType::BlackScholesType
    end


BlackScholesMertonProcess
~~~~~~~~~~~~~~~~~~~~~~~~~

Merton (1973) extension to the Black-Scholes stochastic process

.. function:: BlackScholesMertonProcess(x0::Quote, riskFreeRate::YieldTermStructure, dividendYield::YieldTermStructure, blackVolatility::BlackConstantVol, disc::AbstractDiscretization = EulerDiscretization())

    Constructs a Black Scholes Merton Process, based off the GeneralizedBlackScholesProcess structure above.


Other Stochastic Processes
--------------------------

OrnsteinUhlenbeckProcess
~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type OrnsteinUhlenbeckProcess <: StochasticProcess1D
      speed::Float64
      vol::Float64
      x0::Float64
      level::Float64
    end

.. function:: OrnsteinUhlenbeckProcess(speed::Float64, vol::Float64, x0::Float64 = 0.0, level::Float64 = 0.0) = new(speed, vol, x0, level)

    Constructor for the OrnsteinUhlenbeckProcess


GsrProcess
~~~~~~~~~~

Gaussian short rate stochastic process

.. code-block:: julia

    type GsrProcess <: StochasticProcess1D
      times::Vector{Float64}
      vols::Vector{Float64}
      reversions::Vector{Float64}
      T::Float64
      revZero::Vector{Bool}
      cache1::Dict{Pair{Float64, Float64}, Float64}
      cache2::Dict{Pair{Float64, Float64}, Float64}
      cache3::Dict{Pair{Float64, Float64}, Float64}
      cache4::Dict{Float64, Float64}
      cache5::Dict{Pair{Float64, Float64}, Float64}
    end


.. function:: GsrProcess(times::Vector{Float64}, vols::Vector{Float64}, reversions::Vector{Float64}, T::Float64 = 60.0)

    Constructor for the GsrProcess
