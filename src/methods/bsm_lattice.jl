type BlackScholesLattice{T <: AbstractBinomialTree} <: TreeLattice
  tree::T
  riskFreeRate::Float64
  dt::Float64
  discountFactor::Float64
  pd::Float64
  pu::Float64
  treeLattice::TreeLattice1D

  function BlackScholesLattice(tree::T, riskFreeRate::Float64, endTime::Float64, steps::Int)
    dt = endTime / steps
    discountFactor = exp(-riskFreeRate * dt)
    pd = probability(tree, 0, 0, 0)
    pu = probability(tree, 0, 0, 1)
    bsm_lattice = new{T}(tree, riskFreeRate, dt, discountFactor, pd, pu)

    tg = TimeGrid(endTime, steps)
    lat = TreeLattice1D(tg, 2, bsm_lattice)
    bsm_lattice.treeLattice = lat

    return bsm_lattice
  end
end

get_underlying(bsm::BlackScholesLattice, i::Int, idx::Int) = get_underlying(bsm.tree, i, idx)

descendant(bsm::BlackScholesLattice, i::Int, idx::Int, branch::Int) = descendant(bsm.tree, i, idx, branch)

probability(bsm::BlackScholesLattice, i::Int, idx::Int, branch::Int) = probability(bsm.tree, i, idx, branch)

get_size(bsm::BlackScholesLattice, i::Int) = get_size(bsm.tree, i)

discount(bsm::BlackScholesLattice, ::Int, ::Int) = bsm.discountFactor

function step_back!(bsm::BlackScholesLattice, i::Int, values::Vector{Float64}, newValues::Vector{Float64})
  for j = 1:i
    newValues[j] = (bsm.pd * values[j] + psm.pu * values[j+1]) * bsm.discountFactor
  end

  return newValues
end
