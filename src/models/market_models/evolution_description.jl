type EvolutionDescription
  numberOfRates::Int
  rateTimes::Vector{Float64}
  evolutionTimes::Vector{Float64}
  relevanceRates::Vector{Pair{Int, Int}}
  rateTaus::Vector{Float64}
  firstAliveRate::Vector{Int}
end

function EvolutionDescription(rateTimes::Vector{Float64},
                              evolutionTimes::Vector{Float64},
                              relevanceRates::Vector{Pair{Int, Int}} = Vector{Pair{Int, Int}}())

  numberOfRates = isempty(rateTimes) ? 0 : length(rateTimes) - 1
  evolutionTimes = isempty(evolutionTimes) && ~isempty(rateTimes) ? rateTimes[1:end-1] : evolutionTimes

  firstAliveRate = zeros(length(evolutionTimes))

  rateTaus = check_increasing_times_and_calculate_taus(rateTimes)
  check_increasing_times(evolutionTimes)
  numberOfSteps = length(evolutionTimes)

  evolutionTimes[end] <= rateTimes[end - 1] || error("The last evolution time is past the last fixing time")

  if isempty(relevanceRates)
    relevanceRates = fill(Pair(0, numberOfRates), numberOfSteps)
  else
    length(relevanceRates) == numberOfSteps || error("relevanceRates / evolutionTimes mismatch")
  end

  currentEvolutionTime = 0.0
  firstAliveRateNum = 1
  for j in eachindex(firstAliveRate)
    while rateTimes[firstAliveRateNum] <= currentEvolutionTime
      firstAliveRateNum += 1
    end

    firstAliveRate[j] = firstAliveRateNum
    currentEvolutionTime = evolutionTimes[j]
  end

  return EvolutionDescription(numberOfRates, rateTimes, evolutionTimes, relevanceRates, rateTaus, firstAliveRate)
end
