mutable struct BlackKarasinski{TermStructureConsistentModelType, T <: TermStructure} <: OneFactorModel{TermStructureConsistentModelType}
  modT::TermStructureConsistentModelType
  a::ConstantParameter
  sigma::ConstantParameter
  ts::T
  privateConstraint::PrivateConstraint{ConstantParameter}
  common::ModelCommon
end

function BlackKarasinski(ts::T, a::Float64 = 0.1, sigma = 0.1) where {T <: TermStructure}
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  return BlackKarasinski{TermStructureConsistentModelType, T}(TermStructureConsistentModelType(), a_const, sigma_const, ts, privateConstraint, ModelCommon()) # ShortRateModelCommon())
end

generate_arguments!(m::BlackKarasinski) = m # do nothing

mutable struct BlackKarasinskiDynamics{P <: Parameter} <: ShortRateDynamics
  process::OrnsteinUhlenbeckProcess
  fitting::P
  a::Float64
  sigma::Float64
end

# Dynamics #
BlackKarasinskiDynamics(fitting::P, a::Float64, sigma::Float64) where {P <: Parameter} = BlackKarasinskiDynamics{P}(OrnsteinUhlenbeckProcess(a, sigma), fitting, a, sigma)

short_rate(dynamic::BlackKarasinskiDynamics, t::Float64, x::Float64) = exp(x + dynamic.fitting(t))

struct BlackKarasinskiHelper <: Function
  sz::Int
  xMin::Float64
  dx::Float64
  dt::Float64
  discountBondPrice::Float64
  statePrices::Vector{Float64}
end

function BlackKarasinskiHelper(i::Int, xMin::Float64, dx::Float64, discountBond::Float64, _tree::OneFactorShortRateTree)
  sz = get_size(_tree, i)
  statePrices = get_state_prices!(_tree, i)
  dt = _tree.tg.dt[i]

  return BlackKarasinskiHelper(sz, xMin, dx, dt, discountBond, statePrices)
end

function (bkh::BlackKarasinskiHelper)(theta::Float64)
  val = bkh.discountBondPrice
  x = bkh.xMin

  @simd for j = 1:bkh.sz
    disc = exp(-exp(theta + x) * bkh.dt)
    val -= bkh.statePrices[j] * disc
    x += bkh.dx
  end

  return val
end

# function operator(bkh::BlackKarasinskiHelper)
#   function _inner(theta::Float64)
#     val = bkh.discountBondPrice
#     x = bkh.xMin
#
#     for j = 1:bkh.sz
#       disc = exp(-exp(theta + x) * bkh.dt)
#       val -= bkh.statePrices[j] * disc
#       x += bkh.dx
#     end
#
#     return val
#   end
#
#   return _inner
# end

# tree methods #
function tree(model::BlackKarasinski, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)

  numericDynamics = BlackKarasinskiDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{typeof(numericDynamics), typeof(trinomial)}(trinomial, numericDynamics, grid)

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
    val = solve(slvr, solverHelper, 1e-7, val, vMin, vMax)

    @inbounds set_params!(phi, grid.times[i], val)
  end

  return numericTree
end
