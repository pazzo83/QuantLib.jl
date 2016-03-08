using StatsBase

abstract AbstractStatistics

abstract StatsType
type GaussianStatsType <: StatsType end

type GenericRiskStatistics{S <: StatsType} <: AbstractStatistics
  statsType::S
  samples::Vector{Float64}
  sampleWeights::StatsBase.WeightVec
  samplesMatrix::Matrix{Float64}
  isSorted::Bool
end

gen_RiskStatistics() = GenericRiskStatistics(GaussianStatsType(), Vector{Float64}(), weights(Vector{Float64}()), zeros(0, 2), false)

typealias RiskStatistics GenericRiskStatistics{GaussianStatsType}

function add_sample!(stat::AbstractStatistics, price::Float64, weight::Float64, idx::Int)
  # stat.samplesMatrix[idx, 1] = price
  # stat.samplesMatrix[idx, 2] = weight
  stat.samplesMatrix[idx, :] = [price weight]
  return stat
end

function adding_data!(stat::AbstractStatistics, sz::Int)
  # add to matrix
  stat.samplesMatrix = vcat(stat.samplesMatrix, zeros(sz, 2))
  # add to samples and sample weights
  append!(stat.samples, zeros(sz))
  stat.sampleWeights = weights(vcat(values(stat.sampleWeights), zeros(sz)))

  return stat
end
