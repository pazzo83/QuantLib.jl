branches(::Type{AbstractBinomialTree}) = 2

type BinomialTreeCommon{I <: Integer}
  x0::Float64
  driftPerStep::Float64
  dt::Float64
  columns::I
end

get_x0(t::AbstractBinomialTree) = t.common.x0
get_drift_per_step(t::AbstractBinomialTree) = t.common.driftPerStep
get_dt(t::AbstractBinomialTree) = t.common.dt

function BinomialTreeCommon(process::StochasticProcess1D, endTime::Float64, steps::Int)
  x0 = get_x0(process)
  dt = endTime/steps
  driftPerStep = drift(process, 0.0, x0) * dt
  columns = steps + 1
  BinomialTreeCommon(x0, driftPerStep, dt, columns)
end

type Joshi4 <: BinomialTreeType end
type LeisenReimer <: BinomialTreeType end
type Tian <: BinomialTreeType end

type BinomialTree{T <: BinomialTreeType} <: AbstractBinomialTree
  up::Float64
  down::Float64
  pu::Float64
  pd::Float64
  common::BinomialTreeCommon
  treeType::T
end

typealias Joshi4BinomialTree BinomialTree{Joshi4}
typealias LeisenReimerBinomialTree BinomialTree{LeisenReimer}
typealias TianBinomialTree BinomialTree{Tian}

type JarrowRudd <: EqualProbabilitiesBinomialTreeType end
type AdditiveEQP <: EqualProbabilitiesBinomialTreeType end

type EqualProbabilitiesBinomialTree{T <: EqualProbabilitiesBinomialTreeType} <: AbstractBinomialTree
  up::Float64
  common::BinomialTreeCommon
  treeType::T
end

typealias JarrowRuddBinomialTree EqualProbabilitiesBinomialTree{JarrowRudd}
typealias AdditiveEQPBinomialTree EqualProbabilitiesBinomialTree{AdditiveEQP}

# overloaded to build this alias
function JarrowRudd(process::StochasticProcess1D, endTime::Float64, steps::Int, ::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  up = std_deviation(process, 0.0, common.x0, common.dt)

  return EqualProbabilitiesBinomialTree{JarrowRudd}(up, common, JarrowRudd())
end

type CoxRossRubinstein <: EqualJumpsBinomialTreeType end
type Trigeorgis <: EqualJumpsBinomialTreeType end

type EqualJumpsBinomialTree{T <: EqualJumpsBinomialTreeType} <: AbstractBinomialTree
  dx::Float64
  pu::Float64
  pd::Float64
  common::BinomialTreeCommon
  treeType::T
end

typealias CoxRossRubinsteinBinomialTree EqualJumpsBinomialTree{CoxRossRubinstein}
typealias TrigeorgisBinomialTree EqualJumpsBinomialTree{Trigeorgis}

get_size(tree::AbstractBinomialTree, idx::Int) = idx

descendant(tree::AbstractBinomialTree, ::Int, idx::Int, branch::Int) = (idx - 1) + (branch - 1) + 1

function get_underlying(tree::EqualProbabilitiesBinomialTree, i::Int, idx::Int)
  j = 2 * (idx - 1) - (i - 1)
  # exploiting the forward value tree centering
  return get_x0(tree) * exp((i-1) * get_drift_per_step(tree) + j * tree.up)
end
function get_underlying(tree::EqualJumpsBinomialTree, i::Int, idx::Int)
  j = 2 * idx - i
  # exploting equal jump and the x0 tree centering
  return get_x0(tree) * exp(j * tree.dx)
end
get_underlying(tree::BinomialTree, i::Int, idx::Int) = get_x0(tree) * ^(tree.down, i - idx) * ^(tree.up, idx)

probability(tree::EqualProbabilitiesBinomialTree, ::Int, ::Int, ::Int) = 0.5
probability(tree::Union{EqualJumpsBinomialTree, BinomialTree}, ::Int, ::Int, branch::Int) = branch == 1 ? tree.pu : tree.pd
