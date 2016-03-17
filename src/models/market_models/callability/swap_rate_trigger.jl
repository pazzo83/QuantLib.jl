type SwapRateTrigger <: ExerciseStrategy
  rateTimes::Vector{Float64}
  swapTriggers::Vector{Float64}
  exerciseTimes::Vector{Float64}
  currentIndex::Int
  rateIndex::Vector{Int}
end

function SwapRateTrigger(rateTimes::Vector{Float64}, swapTriggers::Vector{Float64}, exerciseTimes::Vector{Float64})
  check_increasing_times(rateTimes)
  length(rateTimes) > 1 || error("Rate times must contain at least two values")

  check_increasing_times(exerciseTimes)
  length(swapTriggers) == length(exerciseTimes) || error("swapTriggers / exerciseTimes mismatch")

  rateIndex = Vector{Int}(length(exerciseTimes))
  j = 1
  for i in eachindex(exerciseTimes)
    while j <= length(rateTimes) && rateTimes[j] < exerciseTimes[i]
      j += 1
    end

    rateIndex[i] = j
  end

  return SwapRateTrigger(rateTimes, swapTriggers, exerciseTimes, -1, rateIndex)
end

relevant_times(srt::SwapRateTrigger) = srt.exerciseTimes

clone(srt::SwapRateTrigger) = SwapRateTrigger(copy(srt.rateTimes), copy(srt.swapTriggers), copy(srt.exerciseTimes), srt.currentIndex, copy(srt.rateIndex))
