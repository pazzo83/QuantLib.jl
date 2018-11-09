mutable struct SwapForwardBasisSystem <: MarketModelBasisSystem
  rateTimes::Vector{Float64}
  exerciseTimes::Vector{Float64}
  currentIndex::Int
  rateIndex::Vector{Int}
  evolution::EvolutionDescription
end

function SwapForwardBasisSystem(rateTimes::Vector{Float64}, exerciseTimes::Vector{Float64})
  evolution = EvolutionDescription(rateTimes, exerciseTimes)
  rateIndex = Vector{Int}(undef, length(exerciseTimes))
  j = 1
  @inbounds @simd for i in eachindex(exerciseTimes)
    while j <= length(rateTimes) && rateTimes[j] < exerciseTimes[i]
      j += 1
    end

    rateIndex[i] = j
  end

  return SwapForwardBasisSystem(rateTimes, exerciseTimes, -1, rateIndex, evolution)
end

reset!(bs::SwapForwardBasisSystem) = bs.currentIndex = 1

next_step!(bs::SwapForwardBasisSystem, ::CurveState) = bs.currentIndex += 1

function number_of_functions(bs::SwapForwardBasisSystem)
  n = length(bs.exerciseTimes)
  sizes = fill(10, n)
  if bs.rateIndex[n] == length(bs.rateIndex) - 3
    sizes[end] = 6
  end
  if bs.rateIndex[n] == length(bs.rateIndex) - 2
    sizes[end] = 3
  end

  return sizes
end

function set_values!(bs::SwapForwardBasisSystem, currentState::CurveState, results::Vector{Float64})
  rateIndex = bs.rateIndex[bs.currentIndex - 1]

  if rateIndex < length(bs.rateTimes) - 2
    resize!(results, 10)

    x = forward_rate(currentState, rateIndex)
    y = coterminal_swap_rate(currentState, rateIndex + 1)
    z = discount_ratio(currentState, rateIndex, length(bs.rateTimes))

    results[1] = 1.0
    results[2] = x
    results[3] = y
    results[4] = z
    results[5] = x*y
    results[6] = y*z
    results[7] = z*x
    results[8] = x*x
    results[9] = y*y
    results[10] = z*z
  else
    if rateIndex == length(bs.rateTimes) - 2
      x = forward_rate(currentState, rateIndex)
      y = forward_rate(currentState, rateIndex + 1)

      resize!(results, 6)

      results[1] = 1.0
      results[2] = x
      results[3] = y
      results[4] = x*x
      results[5] = x*y
      results[6] = y*y
    else
      x = forward_rate(currentState, rateIndex)

      resize!(results, 3)
      results[1] = 1.0
      results[2] = x
      results[3] = x*x
    end
  end

  return results
end

clone(bs::SwapForwardBasisSystem) = SwapForwardBasisSystem(copy(bs.rateTimes), copy(bs.exerciseTimes), bs.currentIndex, copy(bs.rateIndex), clone(bs.evolution))
