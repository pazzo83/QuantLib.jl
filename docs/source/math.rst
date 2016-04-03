Math Module
===========

This is a sub-module in QuantLib.jl that has math-related types and methods.  It includes various interpolation methods, solvers, and optimization methods.

General Math methods
--------------------

.. function:: is_close{T <: Number}(x::T, y::T, n::Int = 42)

    Determines whether two numbers are almost equal

.. function:: close_enough{T <: Number}(x::T, y::T, n::Int = 42)

    Determines whether two numbers are almost equal but with slightly looser criteria than is_close

.. function:: bounded_log_grid(xMin::Float64, xMax::Float64, steps::Int)

    Bounded log grid constructor

Bernstein  Polynomial
~~~~~~~~~~~~~~~~~~~~~

Bernstein Polynomial type

.. code-block:: julia

    type BernsteinPolynomial end

.. function:: get_polynomial(::BernsteinPolynomial, i::Int, n::Int, x::Float64)

    Build and return a Bernstein polynomial


Interpolation
-------------

QuantLib.jl provides various interpolation methods for building curves

.. code-block:: julia

    abstract Interpolation
    abstract Interpolation2D <: Interpolation

.. function:: update!(interp::Interpolation, idx::Int, val::Float64)

    Updates the y_vals of an interpolation with a given value at a given index


Linear Interpolation
~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type LinearInterpolation <: Interpolation
      x_vals::Vector{Float64}
      y_vals::Vector{Float64}
      s::Vector{Float64}
    end

.. function:: initialize!(interp::LinearInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})

    Initializes the linear interpolation with a given set of x values and y values

.. function:: update!(interp::LinearInterpolation, idx::Int)

    Updates the linear interpolation from a given index

.. function:: update!(interp::LinearInterpolation)

    Updates the linear interpolation from the first index

.. function:: value(interp::LinearInterpolation, val::Float64)

    Returns the interpolated value

.. function:: derivative(interp::LinearInterpolation, val::Float64)

    Returns the derivative of the interpolated value


Log Interpolation
~~~~~~~~~~~~~~~~~

Log interpolation between points - this must be associated with another interpolation method (e.g. Linear).

.. code-block:: julia

    type LogInterpolation <: Interpolation
      x_vals::Vector{Float64}
      y_vals::Vector{Float64}
      interpolator::Interpolation
    end

    typealias LogLinearInterpolation LogInterpolation{LinearInterpolation}

.. function:: LogLinear(x_vals::Vector{Float64}, y_vals::Vector{Float64})

    Constructor for a log linear interpolation

.. function:: LogLinear()

    Construct for a log linear interpolation with no initial values provided

.. function:: initialize!(interp::LogInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})

    Initialize a log interpolation and its interpolator with a given set of x and y values

.. function:: update!(interp::LogInterpolation, idx::Int)

    Updates a log interpolation and its interpolator from a given index

.. function:: update!(interp::LogInterpolation)

    Updates a log interpolation and its interpolator from the first index

.. function:: value(interp::LogInterpolation, val::Float64)

    Returns the interpolated value

.. function:: derivative(interp::LogInterpolation, val::Float64)

    Returns the derivative of the interpolated value


Spline Interpolation
~~~~~~~~~~~~~~~~~~~~

A base cubic interpolation type is provided, with just cubic spline interpolation provided at this time.

.. code-block:: julia
    abstract DerivativeApprox
    abstract BoundaryCondition

    type Spline <: DerivativeApprox end
    type Lagrange <: BoundaryCondition end

    type CubicInterpolation <: Interpolation
      derivativeApprox::DerivativeApprox
      leftBoundaryCondition::BoundaryCondition
      rightBoundaryCondition::BoundaryCondition
      leftValue::Float64
      rightValue::Float64
      monotonic::Bool
      x_vals::Vector{Float64}
      y_vals::Vector{Float64}
      a::Vector{Float64}
      b::Vector{Float64}
      c::Vector{Float64}
      tmp::Vector{Float64}
      dx::Vector{Float64}
      S::Vector{Float64}
      n::Int
      L::TridiagonalOperator
    end

    typealias SplineCubicInterpolation{D, B1, B2} CubicInterpolation{Spline, B1, B2} # First derivative approximation


.. function:: CubicInterpolation(dApprox::DerivativeApprox, leftBoundary::BoundaryCondition, rightBoundary::BoundaryCondition, x_vals::Vector{Float64}, y_vals::Vector{Float64}, leftValue::Float64 = 0.0, rightValue::Float64 = 0.0, monotonic::Bool = true)

    Constructor for any cubic interpolation type

.. function:: value(interp::CubicInterpolation, x::Float64)

    Returns the interpolated value


QuantLib.jl also has a Bicubic spline type, using the Dierckx library.

.. code-block:: julia

    type BicubicSpline
      spline::Dierckx.Spline2D
    end


Backward-Flat Interpolation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

A backward-flat interpolation between points

.. code-block:: julia

    type BackwardFlatInterpolation <: Interpolation
      x_vals::Vector{Float64}
      y_vals::Vector{Float64}
      primitive::Vector{Float64}
    end

.. function:: BackwardFlatInterpolation()

    Constructor for a backward-flat interpolation - no passed in values

.. function:: initialize!(interp::BackwardFlatInterpolation, x_vals::Vector{Float64}, y_vals::Vector{Float64})

    Initializes the backward-flat interpolation with a set of x and y values

.. function:: update!(interp::BackwardFlatInterpolation, idx::Int)

    Updates the backward-flat interpolation from a given index

.. function:: update!(interp::BackwardFlatInterpolation)

    Updates the backward-flat interpolation from the first index

.. function:: value(interp::BackwardFlatInterpolation, x::Float64)

    Returns the interpolated value

.. function:: get_primitive(interp::BackwardFlatInterpolation, x::Float64, ::Bool)

    Returns the primative of the interpolated value


Optimization
------------

QuantLib.jl has various optimization methods for root finding and other calculations

General Optimization types and methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    abstract OptimizationMethod
    abstract CostFunction
    abstract Constraint

OptimizationMethod is the abstract base type.  The CostFunction abstract type is the base type for any function passed to an optimization method.  And the Constraint abstract type provides constraints for the optimization method.

**Projection**

.. code-block:: julia

    type Projection
      actualParameters::Vector{Float64}
      fixedParameters::Vector{Float64}
      fixParams::BitArray
      numberOfFreeParams::Int
    end

The Projection type provides a data structure for actual and fixed parameters for any cost function.

.. function:: project(proj::Projection, params::Vector{Float64})

    Returns the subset of free parameters corresponding to set of parameters

.. function:: include_params(proj::Projection, params::Vector{Float64})

    Returns whole set of parameters corresponding to the set of projected parameters


**Constraints:**

.. code-block:: julia

    type NoConstraint <: Constraint end
    type PositiveConstraint <: Constraint end
    type BoundaryConstraint <: Constraint
      low::Float64
      high::Float64
    end

    type ProjectedConstraint{C <: Constraint} <: Constraint
      constraint::C
      projection::Projection
    end

Various constraint types used in optimization methods


**End Criteria:**

.. code-block:: julia

    type EndCriteria
      maxIterations::Int
      maxStationaryStateIterations::Int
      rootEpsilon::Float64
      functionEpsilon::Float64
      gradientNormEpsilon::Float64
    end

This type provides the end criteria for an optimization method.

**Problem:**

.. code-block:: julia

    type Problem{T}
      costFunction::CostFunction
      constraint::Constraint
      initialValue::Vector{T}
      currentValue::Vector{T}
      functionValue::Float64
      squaredNorm::Float64
      functionEvaluation::Int
      gradientEvaluation::Int
    end

Data structure for a constrained optimization problem.

.. function:: value!{T}(p::Problem, x::Vector{T})

    Calls cost function computation and increment evaluation counter

.. function:: values!(p::Problem, x::Vector{Float64})

    Calls cost values computation and increment evaluation counter


Levenberg Marquardt
~~~~~~~~~~~~~~~~~~~

This is QuantLib.jl's Levenberg Marquardt optimization method

.. code-block:: julia
    type LevenbergMarquardt <: OptimizationMethod
      epsfcn::Float64
      xtol::Float64
      gtol::Float64
      useCostFunctionsJacobin::Bool
    end

.. function:: LevenbergMarquardt()

    Base constructor for the Levenberg Marquardt optimization method

.. function:: minimize!(lm::LevenbergMarquardt, p::Problem, endCriteria::EndCriteria)

    Minimization method for the Levenberg Marquardt optimization method

.. function:: lmdif!(m::Int, n::Int, x::Vector{Float64}, fvec::Vector{Float64}, ftol::Float64, xtol::Float64, gtol::Float64, maxFev::Int, epsfcn::Float64, diag_::Vector{Float64}, mode::Int, factor_::Float64, nprint::Int, info_::Int, nfev::Int, fjac::Matrix{Float64}, ldfjac::Int, ipvt::Vector{Int}, qtf::Vector{Float64}, wa1::Vector{Float64}, wa2::Vector{Float64}, wa3::Vector{Float64}, wa4::Vector{Float64}, fcn!::Function)

    Lmdif function from the MINPACK minimization routine, which has been rewritten in Julia for QuantLib.jl


Simplex
~~~~~~~

This is QuantLib.jl's Simplex optimization method

.. code-block:: julia
    type Simplex <: OptimizationMethod
      lambda::Float64
    end

.. function:: minimize!(simplex::Simplex, p::Problem, end_criteria::EndCriteria)

    Minimization method for the simplex optimization method


Solvers
-------

QuantLib.jl also has several solvers available to find x such that f(x) == 0.


General solver methods
~~~~~~~~~~~~~~~~~~~~~~

These solve methods will call an underlying solve method based on what type of solver you are using.

.. function:: solve(solver::Solver1D, f::Function, accuracy::Float64, guess::Float64, step::Float64)

    General solve method given a guess and step (bounds enforcement is calculated)

.. function:: solve(solver::Solver1D, f::Function, accuracy::Float64, guess::Float64, xMin::Float64, xMax::Float64)

    General solve method given a guess, min, and max.

Mixin shared by all solvers:

.. code-block:: julia

    type SolverInfo
      maxEvals::Int
      lowerBoundEnforced::Bool
      upperBoundEnforced::Bool
      lowerBound::Float64
      upperBound::Float64
    end


Brent Solver
~~~~~~~~~~~~

.. code-block:: julia

    type BrentSolver <: Solver1D
      solverInfo::SolverInfo
    end

.. function:: BrentSolver(maxEvals::Int = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0)

    Constructor for a brent solver, with defaults


Finite Differences Solver
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type FiniteDifferenceNewtonSafe <: Solver1D
      solverInfo::SolverInfo
    end

.. function:: FiniteDifferenceNewtonSafe(maxEvals::Int = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0)

    Constructor for a finite differences newton safe solver, with defaults


Newton Solver
~~~~~~~~~~~~~

.. code-block:: julia

    type NewtonSolver <: Solver1D
      solverInfo::SolverInfo
    end

.. function:: NewtonSolver(maxEvals::Int = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0)

    Constructor for a newton solver, with defaults


General Linear Least Squares
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type GeneralLinearLeastSquares
      a::Vector{Float64}
      err::Vector{Float64}
      residuals::Vector{Float64}
      standardErrors::Vector{Float64}
    end

.. function:: GeneralLinearLeastSquares{T}(x::Vector{Float64}, y::Vector{Float64}, v::Vector{T})

    Constructor for the general linear least squares solver

.. function:: get_coefficients(glls::GeneralLinearLeastSquares) = glls.a

    Returns the coefficients from the general linear least squares calculation


Distributions
-------------

QuantLib.jl largely uses the StatsFuns.jl and Distributions.jl packages for distribution functionality, but we have added a couple additional methods that are necessary for pricing.

.. function:: distribution_derivative(w::Normal, x::Float64)

    Returns the distribution derivative from a normal distribution

.. function:: peizer_pratt_method_2_inversion(z::Float64, n::Int)

    Given an odd integer n and real number z, returns p such that: 1 - CumulativeBinomialDistribution((n-1) / 2, n, p) = CumulativeNormalDistribution(z)


Integration
-----------

QuantLib.jl has various integration methods available.  All are derived from an abstract type Integrator

.. code-block:: julia
    abstract Integrator
    abstract IntegrationFunction

Any integrator has a base "call" method (in Julia, you can make types callable), defined as follows:

.. code-block:: julia

    function call(integrator::Integrator, f::IntegrationFunction, a::Float64, b::Float64)
      integrator.evals = 0
      if a == b
        return 0.0
      end

      if (b > a)
        return integrate(integrator, f, a, b)
      else
        return -integrate(integrator, f, b, a)
      end
    end

Gauss Laguerre Integration
~~~~~~~~~~~~~~~~~~~~~~~~~~

Performs a one-dimensional Gauss-Laguerre integration.

.. code-block:: julia

    abstract GaussianQuadrature

    type GaussLaguerreIntegration <: GaussianQuadrature
      x::Vector{Float64}
      w::Vector{Float64}
    end

The Guass Laguerre integration has its own call method:

.. function:: call(gli::GaussLaguerreIntegration, f::IntegrationFunction)

    This method is called like this, given an instance of the GaussLaguerreIntegration myGLI: myGLI(f) where f is your IntegrationFunction

.. function:: GaussLaguerreIntegration(n::Int, s::Float64 = 0.0)

    Constructor for Gauss Laguerre Integration

.. function:: get_order(integration::GaussianQuadrature)

    Returns the order of the gaussian quadrature


Segment Interval
~~~~~~~~~~~~~~~~

.. code-block:: julia

    type SegmentIntegral <: Integrator
      absoluteAccuracy::Float64
      absoluteError::Float64
      maxEvals::Int
      evals::Int
      intervals::Int
    end

.. function:: SegmentIntegral(intervals::Int)

    Constructor for the segment tntegral

.. function:: Math.integrate(integrator::SegmentIntegral, f::Union{Function, IntegrationFunction}, a::Float64, b::Float64)

    Integrator method for the segment integral


Matrices
--------

Julia has a plethora of matrix-related methods, so we have just added a few specific additional ones needed by various QuantLib.jl calculations

Some specific types used for rank calculation

.. code-block:: julia

    abstract SalvagingAlgo
    type NoneSalvagingAlgo <: SalvagingAlgo end

.. function:: rank_reduced_sqrt(matrix::Matrix, maxRank::Int, componentRetainedPercentage::Float64, sa::NoneSalvagingAlgo)

    Returns the rank-reduced pseudo square root of a real symmetric matrix.  The result matrix has rank<=maxRank. If maxRank>=size, then the specified percentage of eigenvalues out of the eigenvalues' sum is retained.
    If the input matrix is not positive semi definite, it can return an approximation of the pseudo square root using a (user selected) salvaging algorithm.
    The given matrix must be symmetric.


Orthogonal Projection
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: julia

    type OrthogonalProjection
      originalVectors::Matrix{Float64}
      multiplierCutoff::Float64
      numberVectors::Int
      numberValidVectors::Int
      dimension::Int
      validVectors::BitArray{1}
      projectedVectors::Vector{Vector{Float64}}
      orthoNormalizedVectors::Matrix{Float64}
    end

.. function:: OrthogonalProjection(originalVectors::Matrix{Float64}, multiplierCutoff::Float64, tolerance::Float64)

    Given a collection of vectors, w_i, find a collection of vectors x_i such that x_i is orthogonal to w_j for i != j, and <x_i, w_i> = <w_i, w_i>.  This is done by performing GramSchmidt on the other vectors and then projecting onto the orthogonal space.

.. function:: get_vector(op::OrthogonalProjection, i::Int)

    Returns the projected vector of a given index


Random Number Generation
------------------------

We use Julia's built in MersenneTwister RNG for random number generation (RNG), and we have defined types to generate sequences (RSG) based on a few different methods.  We use the StatsFuns and Sobol Julia packages here.

All sequence generators are derived from:

.. code-block:: julia

    abstract AbstractRandomSequenceGenerator

Some general methods:

.. function:: init_sequence_generator!(rsg::AbstractRandomSequenceGenerator, dimension::Int)

    Initializes a RSG with a given dimension (the MersenneTwister, if used, is init-ed with a seed upon type instantiation)

.. function:: last_sequence(rsg::AbstractRandomSequenceGenerator)

    Returns the last sequence generated (the last sequence generated is always cached)


Pseudo Random RSG
~~~~~~~~~~~~~~~~~

This is the most basic random sequence generator.  It simply generates a random sequence based on a set dimension from the Mersenne Twister RNG.

.. code-block:: julia

    type PseudoRandomRSG <: AbstractRandomSequenceGenerator
      rng::MersenneTwister
      dimension::Int
      values::Vector{Float64}
      weight::Float64
    end

.. function:: PseudoRandomRSG(seed::Int, dimension::Int = 1, weight::Float64 = 1.0)

    Constructor for the Pseudo Random RSG, defaults to a sequence length of 1 and weight of 1

.. function:: next_sequence!(rsg::PseudoRandomRSG)

    Builds and returns the next sequence (returns a tuple of the sequence and the weight)


Inverse Cumulative RSG
~~~~~~~~~~~~~~~~~~~~~~

The random numbers in this RSG are generated by the Mersenne Twister and manipulated by the inverse CDF of a normal distribution.

.. code-block:: julia

    type InverseCumulativeRSG <: AbstractRandomSequenceGenerator
      rng::MersenneTwister
      dimension::Int
      values::Vector{Float64}
      weight::Float64
    end

.. function:: InverseCumulativeRSG(seed::Int, dimension::Int = 1, weight::Float64 = 1.0)

    Constructor for the Inverse Cumulative RSG, defaulting to a sequence length of 1 and weight of 1

.. function:: next_sequence!(rsg::InverseCumulativeRSG)

    Builds and returns the next sequence (returns a tuple of the sequence and the weight)


Sobol RSG
~~~~~~~~~

This RSG uses the Julia package Sobol to generate sequences from the Sobol RNG.

.. code-block:: julia

    type SobolRSG <: AbstractRandomSequenceGenerator
      rng::Sobol.SobolSeq
      dimension::Int
      values::Vector{Float64}
      weight::Float64
    end

.. function:: SobolRSG(dimension::Int = 1, weight::Float64 = 1.0)

    Constructor for the Sobol RSG, defaulting to a sequence length of 1 and weight of 1 (no seed required)

.. function:: next_sequence!(rsg::SobolRSG)

    Builds and returns the next sequence (returns a tuple of the sequence and the weight)


Sobol Inverse Cumulative RSG
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This RSG uses the Julia package Sobol to generate Sobol sequences, which are then manipulated via the inverse CDF of a normal distribution.

.. code-block:: julia

    type SobolInverseCumulativeRSG <: AbstractRandomSequenceGenerator
      rng::Sobol.SobolSeq
      dimension::Int
      values::Vector{Float64}
      weight::Float64
    end

.. function:: SobolInverseCumulativeRSG(dimension::Int = 1, weight::Float64 = 1.0)

    Constructor for the Sobol inverse cumulative RSG

.. function:: next_sequence!(rsg::SobolInverseCumulativeRSG)

    Builds and returns the next sequence (returns a tuple of the sequence and the weight)


Sampled Curve
-------------

Contains a sampled curve.  Initially, the structure will contain one indexed curve.

.. code-block:: julia

    type SampledCurve
      grid::Vector{Float64}
      values::Vector{Float64}
    end

.. function:: SampledCurve(gridSize::Int)

    Constructor for a sampled curve with a specified grid size

.. function:: SampledCurve()

    Constructor for a sampled curve with no specified grid size

.. function:: get_size(curve::SampledCurve)

    Returns the size of the sampled curve grid

.. function:: set_log_grid!(curve::SampledCurve, min::Float64, max::Float64)

    Builds a log grid and sets it as the sampled curve's grid

.. function:: set_grid!(curve::SampledCurve, grid::Vector{Float64})

    Sets the sampled curve grid

.. function:: Math.sample!{T}(curve::SampledCurve, f::T)

    Samples from the curve given a passed in function (or something callable)

.. function:: value_at_center(curve::SampledCurve)

    Calculates the value at the center of the sampled curve

.. function:: first_derivative_at_center(curve::SampledCurve)

    Calculates the first derivative at the center of the sampled curve

.. function:: second_derivative_at_center(curve::SampledCurve)

    Calculates the second derivative at the center of the sampled curve


Statistics
----------

These types and methods provide statistical functionality, such as storing data (and weighted data) and performing basic stats calculations (mean, std dev, skewness, etc).  Basic stats functions are provided by StatsBase

.. code-block:: julia

    abstract AbstractStatistics

    abstract StatsType
    type GaussianStatsType <: StatsType end

Basic Stats Methods
~~~~~~~~~~~~~~~~~~~

.. function:: error_estimate(stat::AbstractStatistics)

    Calculates the error estimate of the samples

.. function:: weight_sum(stat::AbstractStatistics)

    Calculates the sum of the weights of all the samples

.. function:: sample_num(stat::AbstractStatistics)

    Returns the number of samples

.. function:: stats_mean(stat::AbstractStatistics)

    Calculates the mean of all the samples

.. function:: stats_std_deviation(stat::AbstractStatistics)

    Calculates the standard deviation of all the samples

.. function:: stats_skewness(stat::AbstractStatistics)

    Calculates the skewness of all the samples

.. function:: stats_kurtosis(stat::AbstractStatistics)

    Calculates the kurtosis of all the samples


Non-weighted Statistics
~~~~~~~~~~~~~~~~~~

The simplest of our stats types, NonWeightedStatistics simply stores non-weighted samples

.. code-block:: julia

    type NonWeightedStatistics <: AbstractStatistics
      samples::Vector{Float64}
      isSorted::Bool
    end

.. function:: NonWeightedStatistics()

    Constructor that initializes an empty NonWeightedStatistics object

.. function:: add_sample!(stat::NonWeightedStatistics, price::Float64)

    Adds a sample


Generic Risk Statistics
~~~~~~~~~~~~~~~~~~~~~~~

This is the basic weighted stats collector

.. code-block:: julia

    type GenericRiskStatistics <: AbstractStatistics
      statsType::StatsType
      samples::Vector{Float64}
      sampleWeights::StatsBase.WeightVec
      samplesMatrix::Matrix{Float64}
      isSorted::Bool
    end

    typealias RiskStatistics GenericRiskStatistics{GaussianStatsType}

.. function:: gen_RiskStatistics(dims::Int = 0)

    Constructs and returns a Risk Statistics object (Generic Risk Statistics with a Gaussian stats type)

.. function:: adding_data!(stat::GenericRiskStatistics, sz::Int)

    This prepares the stats collector to accept new samples.  Because we use a matrix to store the samples and their weights, the size of the matrix has to be preallocated with the number of expected samples.  If you are going to then add additional samples, this must be called again with the number of expected additional samples.

.. function:: add_sample!(stat::GenericRiskStatistics, price::Float64, weight::Float64, idx::Int)

    Adds a new sample with a given weight.  The index is required because we are storing data in a pre-allocated matrix.

.. function:: add_sample!(stat::GenericRiskStatistics, price::Float64, idx::Int)

    Adds a new sample with a default weight of 1.  The index is required because we are storing data in a pre-allocated matrix.


Generic Sequence Statistics
~~~~~~~~~~~~~~~~~~~~~~~~~~~

This data structure stores sequences of AbstractStatistics objects.

.. code-block:: julia

    type GenericSequenceStats <: AbstractStatistics
      dimension::Int
      stats::Vector{AbstractStatistics}
      results::Vector{Float64}
      quadraticSum::Matrix{Float64}
    end

.. function:: GenericSequenceStats(dimension::Int, dim2::Int = dimension)

    Constructor for the generic sequence stats type.  "dimension" is for the number of sequences expected, while "dim2" is for the size of each expected sequence.

.. function:: GenericSequenceStats()

    Constructor for the generic sequence stats type, initializes everything to 0

.. function:: reset!(gss::GenericSequenceStats, dimension::Int, dim2::Int = dimension)

    Resets the generic sequence stats object with the provided dimensions.  "dimension" is for the number of sequences expected, while "dim2" is for the size of each expected sequence.

.. function:: adding_data!(stat::GenericSequenceStats, sz::Int, sz2::Int)

    This is required for adding new data to our sampler.  "sz" is for the number of new sequences added while "sz2" is for the size of those sequences.

.. function:: add_sample!(stat::GenericSequenceStats, vals::Vector, idx::Int, weight::Float64 = 1.0)

    Adds a new sequence of values with a given weight.

.. function:: weight_sum(stat::GenericSequenceStats)

    Calculates the sum of the weights in the first sequence

.. function:: stats_mean(stat::GenericSequenceStats)

    Calculates the mean for each sequence and returns a vector of the means.

.. function:: stats_covariance(stat::GenericSequenceStats)

    Calculates the covariance across all the sample sequences and returns a covariance matrix.


Transformed Grid
----------------

This data structure encapsulates an array of grid points.  It is used primariy in PDE calculations.

.. code-block:: julia

    type TransformedGrid
      grid::Vector{Float64}
      transformedGrid::Vector{Float64}
      dxm::Vector{Float64}
      dxp::Vector{Float64}
      dx::Vector{Float64}
    end

.. function:: TransformedGrid(grid::Vector{Float64}, f::Function)

    Constructor for the transformed grid, given a function to transform the points

.. function:: LogGrid(grid::Vector{Float64})

    Builds a transformed grid based on the log of the grid values


Tridiagonal Operator
--------------------

We are using the custom Tridiagonal Operator from the original QuantLib, rewritten in Julia.

.. code-block:: julia

    type TridiagonalOperator
      diagonal::Vector{Float64}
      lowerDiagonal::Vector{Float64}
      upperDiagonal::Vector{Float64}
      temp::Vector{Float64}
      n::Int
    end

.. function:: TridiagonalOperator(n::Int)

    Constructor for the tridiagonal operator given a dimension

.. function:: TridiagonalOperator()

    Constructor for the tridiagonal operator that defaults to being empty

.. function:: TridiagIdentity(n::Int)

    Builds and returns an identity tridiagonal operator

.. function:: set_first_row!(L::TridiagonalOperator, valB::Float64, valC::Float64)

    Sets the first row of the tridiagonal structure

.. function:: set_mid_row!(L::TridiagonalOperator, i::Int, valA::Float64, valB::Float64, valC::Float64)

    Sets a middle row of the tridiagonal structure, given an index

.. function:: set_last_row!(L::TridiagonalOperator, valA::Float64, valB::Float64)

    Sets the last row of the tridiagonal structure

.. function:: solve_for!(L::TridiagonalOperator, rhs::Vector{Float64}, result::Vector{Float64})

    Solve the linear system for a given right-hand side, with the solution in the result vector.

.. function:: solve_for(L::TridiagonalOperator, rhs::Vector{Float64})

    Solve the linear system, using LU factorization (builds a Julia tridiagonal structure), returns the result vector.
