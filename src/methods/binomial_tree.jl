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

function Tian(process::StochasticProcess1D, endTime::Float64, steps::Int, ::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  q = exp(variance(process, 0.0, common.x0, common.dt))
  r = exp(common.driftPerStep) * sqrt(q)

  up = 0.5 * r * q * (q + 1.0 + sqrt(q * q + 2.0 * q - 3.0))
  down = 0.5 * r * q * (q + 1.0 - sqrt(q * q + 2.0 * q - 3.0))

  pu = (r - down) / (up - down)
  pd = 1.0 - pu

  pu <= 1.0 || error("negative probability")
  pu >= 0.0 || error("negative probability")

  return BinomialTree{Tian}(up, down, pu, pd, common, Tian())
end

function LeisenReimer(process::StochasticProcess1D, endTime::Float64, steps::Int, strike::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  strike > 0.0 || error("strike must be positive")
  oddSteps = isodd(steps) ? steps : steps + 1
  variance_ = variance(process, 0.0, common.x0, endTime)
  ermqdt = exp(common.driftPerStep + 0.5 * variance_ / oddSteps)
  d2 = (log(common.x0 / strike) + common.driftPerStep * oddSteps) / sqrt(variance_)

  pu = peizer_pratt_method_2_inversion(d2, oddSteps)
  pd = 1.0 - pu
  pdash = peizer_pratt_method_2_inversion(d2 + sqrt(variance_), oddSteps)

  up = ermqdt * pdash / pu
  down = (ermqdt - pu * up) / (1.0 - pu)

  return BinomialTree{LeisenReimer}(up, down, pu, pd, common, LeisenReimer())
end

function compute_up_prob(k::Float64, dj::Float64)
  alpha = dj / sqrt(8.0)
  # alpha2 = alpha * alpha
  alpha3 = alpha ^ 3
  alpha5 = alpha ^ 5
  alpha7 = alpha ^ 7
  beta = -0.375 * alpha - alpha3
  gamma = (5.0 / 6.0) * alpha5 + (13.0/12.0) * alpha3 + (25.0 / 128.0) * alpha
  delta = -0.1025 * alpha - 0.9285 * alpha3 - 1.43 * alpha5 - 0.5 * alpha7
  p = 0.5
  rootk = sqrt(k)
  p += alpha / rootk
  p += beta / (k * rootk)
  p += gamma / (k * k * rootk)
  # delete next line to get results for j-three tree
  p += delta / (k * k * k * rootk)
  return p
end

function Joshi4(process::StochasticProcess1D, endTime::Float64, steps::Int, strike::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  strike > 0.0 || error("strike must be positive")
  oddSteps = isodd(steps) ? steps : steps + 1
  variance_ = variance(process, 0.0, common.x0, endTime)
  ermqdt = exp(common.driftPerStep + 0.5 * variance_ / oddSteps)
  d2 = (log(common.x0 / strike) + common.driftPerStep * oddSteps) / sqrt(variance_)

  pu = compute_up_prob((oddSteps - 1.0) / 2.0, d2)
  pd = 1.0 - pu
  pdash = compute_up_prob((oddSteps - 1.0) / 2.0, d2 + sqrt(variance_))

  up = ermqdt * pdash / pu
  down = (ermqdt - pu * up) / (1.0 - pu)

  return BinomialTree{Joshi4}(up, down, pu, pd, common, Joshi4())
end

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

function AdditiveEQP(process::StochasticProcess1D, endTime::Float64, steps::Int, ::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  up = -0.5 * common.driftPerStep + 0.5 * sqrt(4.0 * variance(process, 0.0, common.x0, common.dt) - 3.0 * common.driftPerStep * common.driftPerStep)

  return EqualProbabilitiesBinomialTree{AdditiveEQP}(up, common, AdditiveEQP())
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

function CoxRossRubinstein(process::StochasticProcess1D, endTime::Float64, steps::Int, ::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  dx = std_deviation(process, 0.0, common.x0, common.dt)
  pu = 0.5 + 0.5 * common.driftPerStep / dx
  pd = 1.0 - pu

  pu <= 1.0 || error("negative probability")
  pu >= 0.0 || error("negative probability")

  return EqualJumpsBinomialTree{CoxRossRubinstein}(dx, pu, pd, common, CoxRossRubinstein())
end

function Trigeorgis(process::StochasticProcess1D, endTime::Float64, steps::Int, ::Float64)
  common = BinomialTreeCommon(process, endTime, steps)
  dx = sqrt(variance(process, 0.0, common.x0, common.dt) + common.driftPerStep * common.driftPerStep)
  pu = 0.5 + 0.5 * common.driftPerStep / dx
  pd = 1.0 - pu

  pu <= 1.0 || error("negative probability")
  pu >= 0.0 || error("negative probability")

  return EqualJumpsBinomialTree{Trigeorgis}(dx, pu, pd, common, Trigeorgis())
end

get_size(tree::AbstractBinomialTree, idx::Int) = idx

descendant(tree::AbstractBinomialTree, ::Int, idx::Int, branch::Int) = (idx - 1) + (branch - 1) + 1

function get_underlying(tree::EqualProbabilitiesBinomialTree, i::Int, idx::Int)
  j = 2 * (idx - 1) - (i - 1)
  # exploiting the forward value tree centering
  return get_x0(tree) * exp((i-1) * get_drift_per_step(tree) + j * tree.up)
end
function get_underlying(tree::EqualJumpsBinomialTree, i::Int, idx::Int)
  j = 2 * (idx - 1) - (i - 1)
  # exploting equal jump and the x0 tree centering
  return get_x0(tree) * exp(j * tree.dx)
end
get_underlying(tree::BinomialTree, i::Int, idx::Int) = get_x0(tree) * ^(tree.down, (i - 1) - (idx - 1)) * ^(tree.up, (idx - 1))

probability(tree::EqualProbabilitiesBinomialTree, ::Int, ::Int, ::Int) = 0.5
probability(tree::Union{EqualJumpsBinomialTree, BinomialTree}, ::Int, ::Int, branch::Int) = branch == 2 ? tree.pu : tree.pd
