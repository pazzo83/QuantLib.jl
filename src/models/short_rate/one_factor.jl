## ONE FACTOR MODELS ##
type OneFactorShortRateTree{S <: ShortRateDynamics} <: ShortRateTree
  tree::TrinomialTree
  dynamics::S
  tg::TimeGrid
  treeLattice::TreeLattice1D

  function OneFactorShortRateTree{S}(tree::TrinomialTree, dynamics::S, tg::TimeGrid)
    oneFactorTree = new(tree, dynamics, tg)
    oneFactorTree.treeLattice = TreeLattice1D(tg, get_size(tree, 2), oneFactorTree)

    return oneFactorTree
  end
end

get_size{I <: Integer}(tr::OneFactorShortRateTree, i::I) = get_size(tr.tree, i)

function discount(tr::OneFactorShortRateTree, i::Int, idx::Int)
  x = get_underlying(tr.tree, i, idx)
  r = short_rate(tr.dynamics, tr.tg.times[i], x)
  return exp(-r * tr.tg.dt[i])
end

descendant(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = descendant(tr.tree, i, idx, branch)
probability(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = probability(tr.tree, i, idx, branch)

get_params(m::OneFactorModel) = Float64[get_a(m), get_sigma(m)]

type BlackKarasinski{TermStructureConsistentModelType, T <: TermStructure} <: OneFactorModel{TermStructureConsistentModelType}
  modT::TermStructureConsistentModelType
  a::ConstantParameter
  sigma::ConstantParameter
  ts::T
  privateConstraint::PrivateConstraint
  common::ShortRateModelCommon
end

function BlackKarasinski(ts::TermStructure, a::Float64 = 0.1, sigma = 0.1)
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  return BlackKarasinski(TermStructureConsistentModelType(), a_const, sigma_const, ts, privateConstraint, ShortRateModelCommon())
end

generate_arguments!(m::BlackKarasinski) = m # do nothing


type HullWhite{AffineModelType, T <: TermStructure} <: OneFactorModel{AffineModelType}
  modT::AffineModelType
  r0::Float64
  a::ConstantParameter
  sigma::ConstantParameter
  phi::HullWhiteFittingParameter
  ts::T
  privateConstraint::PrivateConstraint
  common::ShortRateModelCommon
end

function HullWhite{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01)
  _rate = forward_rate(ts, 0.0, 0.0, ContinuousCompounding(), NoFrequency())
  r0 = _rate.rate
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  phi  = HullWhiteFittingParameter(a, sigma, ts)

  return HullWhite(AffineModelType(), r0, a_const, sigma_const, phi, ts, privateConstraint, ShortRateModelCommon())
end

## Dynamics ##

type HullWhiteDynamics{P <: Parameter} <: ShortRateDynamics
  process::OrnsteinUhlenbeckProcess
  fitting::P
  a::Float64
  sigma::Float64
end

HullWhiteDynamics{P <: Parameter}(fitting::P, a::Float64, sigma::Float64) = HullWhiteDynamics(OrnsteinUhlenbeckProcess(a, sigma), fitting, a, sigma)

short_rate(dynamic::HullWhiteDynamics, t::Float64, x::Float64) = x + dynamic.fitting(t)

generate_arguments!(m::HullWhite) = m.phi = HullWhiteFittingParameter(get_a(m), get_sigma(m), m.ts)

get_dynamics(m::HullWhite) = HullWhiteDynamics(m.phi, get_a(m), get_sigma(m))

type BlackKarasinskiDynamics{P <: Parameter} <: ShortRateDynamics
  process::OrnsteinUhlenbeckProcess
  fitting::P
  a::Float64
  sigma::Float64
end

BlackKarasinskiDynamics{P <: Parameter}(fitting::P, a::Float64, sigma::Float64) = BlackKarasinskiDynamics(OrnsteinUhlenbeckProcess(a, sigma), fitting, a, sigma)

short_rate(dynamic::BlackKarasinskiDynamics, t::Float64, x::Float64) = exp(x + dynamic.fitting(t))

immutable BlackKarasinskiHelper{I <: Integer}
  sz::I
  xMin::Float64
  dx::Float64
  dt::Float64
  discountBondPrice::Float64
  statePrices::Vector{Float64}
end

function BlackKarasinskiHelper{I <: Integer}(i::I, xMin::Float64, dx::Float64, discountBond::Float64, _tree::OneFactorShortRateTree)
  sz = get_size(_tree, i)
  statePrices = get_state_prices!(_tree, i)
  dt = _tree.tg.dt[i]

  return BlackKarasinskiHelper(sz, xMin, dx, dt, discountBond, statePrices)
end

function operator(bkh::BlackKarasinskiHelper)
  function _inner(theta::Float64)
    val = bkh.discountBondPrice
    x = bkh.xMin

    for j = 1:bkh.sz
      disc = exp(-exp(theta + x) * bkh.dt)
      val -= bkh.statePrices[j] * disc
      x += bkh.dx
    end

    return val
  end

  return _inner
end

## Tree methods ##

function tree(model::HullWhite, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)
  numericDynamics = HullWhiteDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{HullWhiteDynamics}(trinomial, numericDynamics, grid)

  reset_param_impl!(phi)

  @simd for i = 1:length(grid.times) - 1
    @inbounds discountBond = discount(model.ts, grid.times[i + 1])
    statePrices = get_state_prices!(numericTree, i)
    sz = get_size(numericTree, i)
    @inbounds dt = grid.dt[i]
    @inbounds dx = trinomial.dx[i]
    x = get_underlying(trinomial, i, 1)
    val = 0.0
    for j = 1:sz
      @inbounds val += statePrices[j] * exp(-x * dt)
      x += dx
    end
    val = log(val / discountBond) / dt
    @inbounds set_params!(phi, grid.times[i], val)
  end

  return numericTree
end

function tree(model::BlackKarasinski, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)

  numericDynamics = BlackKarasinskiDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{BlackKarasinskiDynamics}(trinomial, numericDynamics, grid)

  reset_param_impl!(phi)

  val = 1.0
  vMin = -50.0
  vMax = 50.0

  @simd for i = 1:length(grid.times) - 1
    @inbounds discountBond = discount(model.ts, grid.times[i + 1])
    xMin = get_underlying(trinomial, i, 1)
    @inbounds dx = trinomial.dx[i]

    solverHelper = BlackKarasinskiHelper(i, xMin, dx, discountBond, numericTree)
    slvr = BrentSolver(1000) # max evals = 1000
    val = solve(slvr, operator(solverHelper), 1e-7, val, vMin, vMax)

    @inbounds set_params!(phi, grid.times[i], val)
  end

  return numericTree
end


type RStarFinder{M <: ShortRateModel}
  model::M
  strike::Float64
  maturity::Float64
  valueTime::Float64
  fixedPayTimes::Vector{Float64}
  amounts::Vector{Float64}
end

function operator(rsf::RStarFinder)
  function _inner(x::Float64)
    _value = rsf.strike
    _B = discount_bond(rsf.model, rsf.maturity, rsf.valueTime, x)
    sz = length(rsf.fixedPayTimes)
    for i = 1:sz
      dbVal = discount_bond(rsf.model, rsf.maturity, rsf.fixedPayTimes[i], x) / _B
      _value -= rsf.amounts[i] * dbVal
    end

    return _value
  end

  return _inner
end

function B(model::HullWhite, t::Float64, T::Float64)
  _a = get_a(model)
  if _a < sqrt(eps())
    return T - t
  else
    return (1.0 - exp(-_a * (T - t))) / _a
  end
end

function A(model::HullWhite, t::Float64, T::Float64)
  discount1 = discount(model.ts, t)
  discount2 = discount(model.ts, T)

  forward = forward_rate(model.ts, t, t, ContinuousCompounding(), NoFrequency())
  temp = get_sigma(model) * B(model, t, T)
  val = B(model, t, T) * forward.rate - 0.25 * temp * temp * B(model, 0.0, 2.0* t)

  return exp(val) * discount2 / discount1
end

# type HullWhite{T <: TermStructure} <: ShortRateModel
#   r0::Float64
#   a::ConstantParameter
#   sigma::ConstantParameter
#   phi::HullWhiteFittingParameter
# end

discount_bond(model::OneFactorModel, tNow::Float64, maturity::Float64, factors::Vector{Float64}) = discount_bond(model, tNow, maturity, factors[1])
discount_bond{M <: OneFactorModel}(model::M, tNow::Float64, maturity::Float64, _rate::Float64) = A(model, tNow, maturity) * exp(-B(model, tNow, maturity) * _rate)

function discount_bond_option{O <: OptionType}(model::HullWhite, optionType::O, strike::Float64, maturity::Float64, bondStart::Float64, bondMaturity::Float64)
  _a = get_a(model)
  if _a < sqrt(eps())
    v = get_sigma(model) * B(model, bondStart, bondMaturity) * sqrt(maturity)
  else
    v = get_sigma(model) / (_a * sqrt(2.0 * _a)) * sqrt(exp(-2.0 * _a * (bondStart - maturity)) - exp(-2.0 * _a * bondStart) - 2.0 * (exp(-_a * (bondStart + bondMaturity - 2.0 * maturity))
        - exp(-_a * (bondStart + bondMaturity))) + exp(-2.0 * _a * (bondMaturity - maturity)) - exp(-2.0 * _a * bondMaturity))
  end

   f = discount(model.ts, bondMaturity)
   k = discount(model.ts, bondStart) * strike

   return black_formula(optionType, k, f, v)
 end
