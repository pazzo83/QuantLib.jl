type MultiStepInverseFloater
  rateTimes::Vector{Float64}
  evolution::EvolutionDescription
  fixedAccruals::Vector{Float64}
  floatingAccruals::Vector{Float64}
  fixedStrikes::Vector{Float64}
  fixedMultipliers::Vector{Float64}
  floatingSpreads::Vector{Float64}
  paymentTimes::Vector{Float64}
  payer::Bool
  multipler::Float64
  lastIndex::Int
  currentIndex::Int
end

function MultiStepInverseFloater(rateTimes::Vector{Float64},
                                fixedAccruals::Vector{Float64},
                                floatingAccruals::Vector{Float64},
                                fixedStrikes::Vector{Float64},
                                fixedMultipliers::Vector{Float64},
                                floatingSpreads::Vector{Float64},
                                paymentTimes::Vector{Float64},
                                payer::Bool = true)
  n = length(rateTimes) - 1
  evolutionTimes = zeros(n)
  relevanceRates = Vector{Pair{Int, Int}}(n)
  for i in eachindex(evolutionTimes)
    evolutionTimes[i] = rateTimes[i]
    relevanceRates[i] = Pair(i, i+1)
  end

  evolution = EvolutionDescription(rateTimes, evolutionTimes, relevanceRates)
  multipler = payer ? -1.0 : 1.0

  return MultiStepInverseFloater(rateTimes, evolution, fixedAccruals, floatingAccruals, fixedStrikes, fixedMultipliers, floatingSpreads, paymentTimes,
                                payer, multipler, n+1, -1)
end
