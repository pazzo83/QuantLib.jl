type NothingExerciseValue
  numberOfExercises::Int
  rateTimes::Vector{Float64}
  isExerciseTime::BitArray{1}
  evolution::EvolutionDescription
  currentIndex::Int
  cf::MarketModelCashFlow
end

function NothingExerciseValue(rateTimes::Vector{Float64}, isExerciseTime::BitArray{1} = BitArray{1}())
  currentIndex = 1
  check_increasing_times(rateTimes)
  length(rateTimes) >= 2 || error("Rate times must contain at least two values")

  cf = MarketModelCashFlow(-1, 0.0)
  evolutionTimes = rateTimes[1:end-1]
  evolution = EvolutionDescription(rateTimes, evolutionTimes)

  if isempty(isExerciseTime)
    isExerciseTime = trues(isempty(rateTimes) ? 0 : length(rateTimes) - 1)
  else
    length(isExerciseTime) == isempty(rateTimes) ? 0 : length(rateTimes) - 1 || error("isExerciseTime must have same size as rateTimes minus 1")
  end

  numberOfExercises = 0
  for i in eachindex(isExerciseTime)
    if isExerciseTime[i]
      numberOfExercises += 1
    end
  end

  return NothingExerciseValue(numberOfExercises, rateTimes, isExerciseTime, evolution, currentIndex, cf)
end
