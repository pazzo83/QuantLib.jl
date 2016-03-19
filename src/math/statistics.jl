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

gen_RiskStatistics(dims::Int = 0) = GenericRiskStatistics(GaussianStatsType(), Vector{Float64}(dims), weights(Vector{Float64}(dims)), zeros(dims, 2), false)

typealias RiskStatistics GenericRiskStatistics{GaussianStatsType}

type GenericSequenceStats{S <: AbstractStatistics} <: AbstractStatistics
  dimension::Int
  stats::Vector{S}
  results::Vector{Float64}
  quadraticSum::Matrix{Float64}
end

function GenericSequenceStats(dimension::Int, dim2::Int = dimension)
  # initial init
  gss = GenericSequenceStats(0, Vector{GenericRiskStatistics}(0), Vector{Float64}(0), Matrix{Float64}(0, 0))

  reset!(gss, dimension, dim2)

  return gss
end

function reset!(gss::GenericSequenceStats, dimension::Int, dim2::Int = dimension)
  # re-init
  if dimension > 0
    if dimension == gss.dimension
      for i in eachindex(gss.stats)
        reset!(gss.stats[i])
      end
    else
      gss.dimension = dimension
      gss.stats = GenericRiskStatistics[gen_RiskStatistics(dim2) for i = 1:dimension]
      results = zeros(dimension)
    end
    gss.quadraticSum = zeros(dimension, dimension)
  else
    gss.dimension = dimension
  end

  return gss
end

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

function add_sample!(stat::GenericSequenceStats, vals::Vector, idx::Int, weight::Float64 = 1.0)
  if stat.dimension == 0
    # stat wasn't initialized
    dimension = length(vals)
    reset!(stat, dimension)
  end

  length(vals) == stat.dimension || error("sample size mismatch")

  BLAS.ger!(1.0, vals, vals, stat.quadraticSum)

  for i in eachindex(stat.stats)
    add_sample!(stat.stats[i], vals[i], weight, idx)
  end

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

function stats_mean(stat::GenericSequenceStats)
  return Float64[stats_mean(x) for x in stat.stats]
end

function stats_std_deviation(stat::AbstractStatistics)
  sort_samples!(stat)
  return std(stat.samples)
end

function stats_skewness(stat::AbstractStatistics)
  sort_samples!(stat)
  return skewness(stat.samples)
end

function stats_kurtosis(stat::AbstractStatistics)
  sort_samples!(stat)
  return kurtosis(stat.samples)
end
