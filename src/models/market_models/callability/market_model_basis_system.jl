type SwapForwardBasisSystem <: MarketModelBasisSystem
  rateTimes::Vector{Float64}
  exerciseTimes::Vector{Float64}
  currentIndex::Int
  rateIndex::Vector{Int}
  evolution::EvolutionDescription
end

function SwapForwardBasisSystem(rateTimes::Vector{Float64}, exerciseTimes::Vector{Float64})
  evolution = EvolutionDescription(rateTimes, exerciseTimes)
  rateIndex = Vector{Int}(length(exerciseTimes))
  j = 1
  for i in eachindex(exerciseTimes)
    while j <= length(rateTimes) && rateTimes[j] < exerciseTimes[i]
      j += 1
    end

    rateIndex[i] = j
  end

  return SwapForwardBasisSystem(rateTimes, exerciseTimes, -1, rateIndex, evolution)
end
