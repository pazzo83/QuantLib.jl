using StatsBase

abstract AbstractStatistics

abstract StatsType
type GaussianStatsType <: StatsType end

type NonWeightedStatistics <: AbstractStatistics
  samples::Vector{Float64}
  isSorted::Bool
end

NonWeightedStatistics() = NonWeightedStatistics(Vector{Float64}(), false)

type GenericRiskStatistics{S <: StatsType} <: AbstractStatistics
  statsType::S
  samples::Vector{Float64}
  sampleWeights::StatsBase.WeightVec
  samplesMatrix::Matrix{Float64}
  isSorted::Bool
end

gen_RiskStatistics() = GenericRiskStatistics(GaussianStatsType(), Vector{Float64}(), weights(Vector{Float64}()), zeros(0, 2), false)

typealias RiskStatistics GenericRiskStatistics{GaussianStatsType}

function add_sample!(stat::GenericRiskStatistics, price::Float64, weight::Float64, idx::Int)
  # stat.samplesMatrix[idx, 1] = price
  # stat.samplesMatrix[idx, 2] = weight
  stat.samplesMatrix[idx, :] = [price weight]
  stat.isSorted = false
  return stat
end

function add_sample!(stat::NonWeightedStatistics, price::Float64)
  push!(stat.samples, price)
  stat.isSorted = false
  return stat
end

function adding_data!(stat::GenericRiskStatistics, sz::Int)
  # add to matrix
  stat.samplesMatrix = vcat(stat.samplesMatrix, zeros(sz, 2))
  # add to samples and sample weights
  append!(stat.samples, zeros(sz))
  stat.sampleWeights = weights(vcat(values(stat.sampleWeights), zeros(sz)))
  stat.isSorted = false
  return stat
end

function sort_samples!(stat::GenericRiskStatistics)
  if ~stat.isSorted
    stat.samplesMatrix = sortrows(stat.samplesMatrix)
    stat.samples = stat.samplesMatrix[:, 1]
    stat.sampleWeights = weights(stat.samplesMatrix[:, 2])
    stat.isSorted = true
  end
  return stat
end

function sort_samples!(stat::NonWeightedStatistics)
  sort!(stat.samples)
  stat.isSorted = true
end

function error_estimate(stat::AbstractStatistics)
  sort_samples!(stat)
  return sqrt(var(stat.samples)/ length(stat.samples))
end

function stats_mean(stat::AbstractStatistics)
  sort_samples!(stat)
  return mean(stat.samples)
end
