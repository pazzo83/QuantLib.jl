Pricing Methods
===============

QuantLib.jl has various methods for asset pricing and calculation.

Finite Differences
------------------

Finite Differences framework for option pricing

General Finite Differences types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

FDM Solver Description:

.. code-block:: julia
    type FdmSolverDesc{F <: FdmMesher, C <: FdmInnerValueCalculator}
      mesher::F
      bcSet::FdmBoundaryConditionSet
      condition::FdmStepConditionComposite
      calculator::C
      maturity::Float64
      timeSteps::Int
      dampingSteps::Int
    end


FDM Calculation Related Types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Inner value calculators

.. code-block:: julia

    abstract FdmInnerValueCalculator

FDM Affine Model Swap Inner Value Calculator:

.. code-block:: julia

    type FdmAffineModelSwapInnerValue <: FdmInnerValueCalculator
      disModel::Model
      fwdModel::Model
      swap::VanillaSwap
      exerciseDates::Dict{Float64, Date}
      mesher::FdmMesher
      direction::Int
      disTs::FdmAffineModelTermStructure
      fwdTs::FdmAffineModelTermStructure
    end

.. function:: FdmAffineModelSwapInnerValue(disModel::Model, fwdModel::Model, swap::VanillaSwap, exerciseDates::Dict{Float64, Date}, mesher::FdmMesher, direction::Int)

    Constructor for the FDM Affine Swap Inner Value calculator, given a discount model, a forward model, a swap, and mesher


Finite Difference Step Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    abstract StepCondition


FDM Composite Step Condition

.. code-block:: julia

    type FdmStepConditionComposite{C <: StepCondition} <: StepCondition
      stoppingTimes::Vector{Float64}
      conditions::Vector{C}
    end

.. function:: vanilla_FdmStepConditionComposite(cashFlow::DividendSchedule, exercise::Exercise, mesher::FdmMesher, calculator::FdmInnerValueCalculator, refDate::Date, dc::DayCount)

    Builds and returns an FDM Step Condition composite for vanilla options


Finite Difference meshers
~~~~~~~~~~~~~~~~~~~~~~~~~

Brief mesher for an FDM grid

.. code-block:: julia

    abstract FdmMesher
    abstract Fdm1DMesher


Composite FDM Mesher type

.. code-block:: julia

    type FdmMesherComposite <: FdmMesher
      layout::FdmLinearOpLayout
      meshers::Vector{Fdm1DMesher} # this could change
    end

.. function:: FdmMesherComposite(mesh::Fdm1DMesher)

    Constructor for FDM Mesher composite type, with one mesher

.. function:: FdmMesherComposite{F1D <: Fdm1DMesher}(xmesher::F1D, ymesher::F1D)

    Constructor for FDM Mesher composite type with two meshers of the same type


1D FDM Mesher using a stochastic process:

.. code-block:: julia

    type FdmSimpleProcess1dMesher <: Fdm1DMesher
      size::Int
      process::StochasticProcess1D
      maturity::Float64
      tAvgSteps::Int
      epsilon::Float64
      mandatoryPoint::Float64
      locations::Vector{Float64}
      dplus::Vector{Float64}
      dminus::Vector{Float64}
    end

.. function:: FdmSimpleProcess1dMesher(sz::Int, process::StochasticProcess1D, maturity::Float64, tAvgSteps::Int, _eps::Float64, mandatoryPoint::Float64 = -1.0)

    Constructor for the FDM Simple Process (stochastic process) 1D mesher

Finite Difference operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Linear operators to model a multi-dimensional PDE system

.. code-block:: julia

    abstract FdmLinearOpComposite


FDM G2 operator - linear operator for the FDM G2 solver

.. code-block:: julia

    type FdmG2Op <: FdmLinearOpComposite
      direction1::Int
      direction2::Int
      x::Vector{Float64}
      y::Vector{Float64}
      dxMap::TripleBandLinearOp
      dyMap::TripleBandLinearOp
      corrMap::SecondOrderMixedDerivativeOp
      mapX::TripleBandLinearOp
      mapY::TripleBandLinearOp
      model::G2
    end

.. function:: FdmG2Op(mesher::FdmMesher, model::G2, direction1::Int, direction2::Int)

    Constructor for the FDM G2 operator

Finite Difference 2D Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~

2D Finite Differences solver

.. code-block:: julia

    type Fdm2DimSolver <: LazyObject
      lazyMixin::LazyMixin
      solverDesc::FdmSolverDesc
      schemeDesc::FdmSchemeDesc
      op::FdmLinearOpComposite
      thetaCondition::FdmSnapshotCondition
      conditions::FdmStepConditionComposite
      initialValues::Vector{Float64}
      resultValues::Matrix{Float64}
      x::Vector{Float64}
      y::Vector{Float64}
      interpolation::BicubicSpline
    end


.. function:: Fdm2DimSolver(solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc, op::FdmLinearOpComposite)

    Constructor for the FDM 2D solver



Finite Difference G2 Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finite differences solver with a G2 short rate model

.. code-block:: julia

    type FdmG2Solver <: LazyObject
      lazyMixin::LazyMixin
      model::G2
      solverDesc::FdmSolverDesc
      schemeDesc::FdmSchemeDesc
      solver::Fdm2DimSolver
    end

.. function:: FdmG2Solver(model::G2, solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc)

    Constructor for the FDM G2 Solver


Finite Difference Hull White Solver
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type FdmHullWhiteSolver <: LazyObject
      lazyMixin::LazyMixin
      model::HullWhite
      solverDesc::FdmSolverDesc
      schemeDesc::FdmSchemeDesc
      solver::Fdm1DimSolver
    end

.. function:: FdmHullWhiteSolver(model::HullWhite, solverDesc::FdmSolverDesc, schemeDesc::FdmSchemeDesc)

    Constructor for the FDM Hull White solver


Monte Carlo Simulation
----------------------

QuantLib.jl's Monte Carlo simulation tools include path generation and path pricer types

Monte Carlo Model
~~~~~~~~~~~~~~~~~

This is the general monte carlo model that is used for the simulation

.. code-block:: julia

  abstract AbstractMonteCarloModel

  type MonteCarloModel <: AbstractMonteCarloModel
    pathGenerator::PathGenerator
    pathPricer::AbstractPathPricer
    sampleAccumulator::RiskStatistics
    isAntitheticVariate::Bool
  end

.. function:: add_samples!(mcmodel::MonteCarloModel, samples::Int, idx::Int=1)

    Adds samples generated by the simulation


Path Generator
~~~~~~~~~~~~~~

.. code-block:: julia

    type PathGenerator
      brownianBridge::Bool
      generator::AbstractRandomSequenceGenerator
      dimension::Int
      timeGrid::TimeGrid
      process::StochasticProcess1D
      nextSample::Sample
      temp::Vector{Float64}
      bb::BrownianBridge
    end

.. function:: PathGenerator(process::StochasticProcess, tg::TimeGrid, generator::AbstractRandomSequenceGenerator, brownianBridge::Bool)

    Constructor for the path generator, given a process, time grid, RSG, and brownian bridge flag

.. function:: PathGenerator(process::StochasticProcess, len::Float64, timeSteps::Int, generator::AbstractRandomSequenceGenerator, brownianBridge::Bool)

    Constructor for the path generator, given a process, two variables to build a time grid, a RSG, and brownian bridge flag


Path and Node
~~~~~~~~~~~~~

Path: Basic path data structure

.. code-block:: julia

    type Path
      tg::TimeGrid
      values::Vector{Float64}
    end

.. function:: Path(tg::TimeGrid)

    Constructor for a Path given a time grid

Node: Node data structure

.. code-block:: julia

    type NodeData
      exerciseValue::Float64
      cumulatedCashFlows::Float64
      values::Vector{Float64}
      controlValue::Float64
      isValid::Bool
    end

.. function:: NodeData()

    Constructor for an empty NodeData object


Misc Methods and Types
~~~~~~~~~~~~~~~~~~~~~~

.. function:: generic_longstaff_schwartz_regression!(simulationData::Vector{Vector{NodeData}}, basisCoefficients::Vector{Vector{Float64}})

    Returns the biased estimate obtained while regressing n exercises, n+1 elements in simulationData
    simulationData[1][j] -> cashflows up to first exercise, j-th path
    simulationData[i+1][j] -> i-th exercise, j-th path
    length(basisCoefficients) = n


Lattices and Trees
------------------

Trinomial Tree
~~~~~~~~~~~~~~

This type defines a recombining trinomial tree approximating a 1D stochastic process.

.. code-block:: julia

    type TrinomialTree <: AbstractTree
      process::StochasticProcess
      timeGrid::TimeGrid
      dx::Vector{Float64}
      branchings::Vector{Branching}
      isPositive::Bool
    end

.. function:: TrinomialTree(process::StochasticProcess, timeGrid::TimeGrid, isPositive::Bool = false)

    Constructor for a trinomial tree given a stochastic process


Tree Lattice 1D
~~~~~~~~~~~~~~~

One-dimensional tree-based lattice

.. code-block:: julia

    type TreeLattice1D <: TreeLattice
      tg::TimeGrid
      impl::TreeLattice
      statePrices::Vector{Vector{Float64}}
      n::Int
      statePricesLimit::Int
    end

.. function:: TreeLattice1D(tg::TimeGrid, n::Int, impl::TreeLattice)

    Constructor of a 1D Tree lattice
