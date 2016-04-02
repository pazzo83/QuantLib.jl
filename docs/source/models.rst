Pricing Models
==============

QuantLib.jl has various pricing models for asset pricing.

Equity Models
-------------

Bates Model
~~~~~~~~~~~

Bates stochastic volatility model

.. code-block:: julia

    type BatesModel{CalibratedModelType} <: Model{CalibratedModelType}
      modT::CalibratedModelType
      hestonModel::HestonModel
      nu::ConstantParameter
      delta::ConstantParameter
      lambda::ConstantParameter
      process::BatesProcess
      common::ModelCommon
    end

.. function:: BatesModel(process::BatesProcess)

    Constructor for the Bates model, given a Bates process

Heston Model
~~~~~~~~~~~~

Heston model for the stochastic volatility of an asset

.. code-block:: julia

    type HestonModel{CalibratedModelType} <: Model{CalibratedModelType}
      modT::CalibratedModelType
      theta::ConstantParameter
      kappa::ConstantParameter
      sigma::ConstantParameter
      rho::ConstantParameter
      v0::ConstantParameter
      process::HestonProcess
      common::ModelCommon
    end

.. function:: HestonModel(process::HestonProcess)

    Constructor for a Heston model given a Heston process


Market Models
-------------
