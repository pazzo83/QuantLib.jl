abstract AbstractStatistics

abstract RiskType
type GaussianRiskType <: RiskType end

type GenericRiskStatistics
  samples::Vector{Float64}
  isSorted::Bool
end
