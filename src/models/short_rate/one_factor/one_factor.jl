## ONE FACTOR MODELS ##
type OneFactorShortRateTree{S <: ShortRateDynamics, P <: StochasticProcess} <: ShortRateTree
  tree::TrinomialTree{P}
  dynamics::S
  tg::TimeGrid
  treeLattice::TreeLattice1D{OneFactorShortRateTree{S, P}}

  function OneFactorShortRateTree{S, P}(tree::TrinomialTree{P}, dynamics::S, tg::TimeGrid)
    oneFactorTree = new{S, P}(tree, dynamics, tg)
    oneFactorTree.treeLattice = TreeLattice1D(tg, get_size(tree, 2), oneFactorTree)

    return oneFactorTree
  end
end

get_size(tr::OneFactorShortRateTree, i::Int) = get_size(tr.tree, i)

function discount(tr::OneFactorShortRateTree, i::Int, idx::Int)
  x = get_underlying(tr.tree, i, idx)
  r = short_rate(tr.dynamics, tr.tg.times[i], x)
  return exp(-r * tr.tg.dt[i])
end

descendant(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = descendant(tr.tree, i, idx, branch)
probability(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = probability(tr.tree, i, idx, branch)

get_params(m::OneFactorModel) = Float64[get_a(m), get_sigma(m)]

type RStarFinder{M <: ShortRateModel} <: Function
  model::M
  strike::Float64
  maturity::Float64
  valueTime::Float64
  fixedPayTimes::Vector{Float64}
  amounts::Vector{Float64}
end

function (rsf::RStarFinder)(x::Float64)
  _value = rsf.strike
  _B = discount_bond(rsf.model, rsf.maturity, rsf.valueTime, x)
  sz = length(rsf.fixedPayTimes)
  @simd for i = 1:sz
    @inbounds dbVal = discount_bond(rsf.model, rsf.maturity, rsf.fixedPayTimes[i], x) / _B
    @inbounds _value -= rsf.amounts[i] * dbVal
  end

  return _value
end

# function operator(rsf::RStarFinder)
#   function _inner(x::Float64)
#     _value = rsf.strike
#     _B = discount_bond(rsf.model, rsf.maturity, rsf.valueTime, x)
#     sz = length(rsf.fixedPayTimes)
#     for i = 1:sz
#       dbVal = discount_bond(rsf.model, rsf.maturity, rsf.fixedPayTimes[i], x) / _B
#       _value -= rsf.amounts[i] * dbVal
#     end
#
#     return _value
#   end
#
#   return _inner
# end

discount_bond(model::OneFactorModel, tNow::Float64, maturity::Float64, factors::Vector{Float64}) = discount_bond(model, tNow, maturity, factors[1])
discount_bond{M <: OneFactorModel}(model::M, tNow::Float64, maturity::Float64, _rate::Float64) = A(model, tNow, maturity) * exp(-B(model, tNow, maturity) * _rate)
