type EvolutionDescription
  numberOfRates::Int
  rateTimes::Vector{Float64}
  evolutionTimes::Vector{Float64}
  relevanceRates::Vector{Pair{Int, Int}}
  rateTaus::Vector{Float64}
  firstAliveRate::Vector{Int}
end

# Constructors #
EvolutionDescription() = EvolutionDescription(-1, Vector{Float64}(), Vector{Float64}(), Vector{Pair{Int, Int}}(), Vector{Float64}(), Vector{Int}())

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

number_of_steps(evol::EvolutionDescription) = length(evol.evolutionTimes)

clone(ed::EvolutionDescription) = EvolutionDescription(ed.numberOfRates, copy(ed.rateTimes), copy(ed.evolutionTimes), copy(ed.relevanceRates),
                                  copy(ed.rateTaus), copy(ed.firstAliveRate))

## Misc functions ##
function money_market_plus_measure(ev::EvolutionDescription, offset::Int)
  rateTimes = ev.rateTimes
  maxNumeraire = length(rateTimes)

  offset <= maxNumeraire || error("offset is greater than the max allowed value for numeraire")
  evolutionTimes = ev.evolutionTimes
  n = length(evolutionTimes)
  numeraires = Vector{Int}(n)
  j = 1
  for i = 1:n
    while rateTimes[j] < evolutionTimes[i]
      j += 1
    end
    numeraires[i] = min(j + offset, maxNumeraire)
  end

  return numeraires
end

money_market_measure(evol::EvolutionDescription) = money_market_plus_measure(evol, 0)

function check_compatibility(evol::EvolutionDescription, numeraires::Vector{Int})
  evolutionTimes = evol.evolutionTimes
  n = length(evolutionTimes)

  length(numeraires) == n || error("size mismatch between numeraires and evolution times")

  rateTimes = evol.rateTimes
  for i = 1:n-1
    rateTimes[numeraires[i]] >= evolutionTimes[i] || error("$(i+1) step, evolution time $(evolutionTimes[i]): the numeraire ($(numeraires[i]))
                                                        corresponding to the rate time $(rateTimes[numeraires[i]]) is expired")
  end
end
