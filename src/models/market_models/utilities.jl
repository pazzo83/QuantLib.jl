type MarketModelCashFlow
  timeIndex::Int
  amount::Float64
end

MarketModelCashFlow() = MarketModelCashFlow(1, 0.0)
clone(mmcf::MarketModelCashFlow) = MarketModelCashFlow(mmcf.timeIndex, mmcf.amount)

type MarketModelPathWiseCashFlow
  timeIndex::Int
  amount::Vector{Float64}
end

MarketModelPathWiseCashFlow(n::Int) = MarketModelPathWiseCashFlow(1, zeros(n))

## Discounter ##
type MarketModelDiscounter
  beforeSize::Int
  beforeWeight::Float64
end

function MarketModelDiscounter(paymentTime::Float64, rateTimes::Vector{Float64})
  check_increasing_times(rateTimes)
  beforeTime = findfirst(rateTimes, paymentTime)

  # handle the case where the payment is in the last period or after the last period
  if beforeTime > length(rateTimes) - 1
    beforeTime = length(rateTimes) - 1
  end

  if beforeTime == 0
    if paymentTime > rateTimes[end]
      beforeTime = length(rateTimes) - 1
    else
      beforeTime = 1
    end
  end

  beforeWeight = 1.0 - (paymentTime - rateTimes[beforeTime]) / (rateTimes[beforeTime + 1] - rateTimes[beforeTime])

  return MarketModelDiscounter(beforeTime, beforeWeight)
end

function numeraire_bonds(mmd::MarketModelDiscounter, curveState::CurveState, numeraire::Int)
  preDF = discount_ratio(curveState, mmd.beforeSize, numeraire)
  if mmd.beforeWeight == 1.0
    return preDF
  end

  postDF = discount_ratio(curveState, mmd.beforeSize+1, numeraire)
  if mmd.beforeWeight == 0.0
    return postDF
  end

  return ^(preDF, mmd.beforeWeight) * ^(postDF, 1.0 - mmd.beforeWeight)
end

function check_increasing_times(times::Vector{Float64})
  length(times) > 0 || error("at least one time is required")
  times[1] > 0.0 || error("first time must be greater than zero")

  issorted(times) || error("non increasing rate times")
end

function check_increasing_times_and_calculate_taus(times::Vector{Float64})
  check_increasing_times(times)

  nTimes = length(times)
  nTimes > 1 || error("at least two times are required")

  taus = Vector{Float64}(nTimes - 1)

  for i in eachindex(taus)
    taus[i] = times[i+1] - times[i]
  end

  return taus
end

function merge_times(times::Vector{Vector{Float64}})
  allTimes = Vector{Float64}()
  for i in eachindex(times)
    allTimes = vcat(allTimes, times[i])
  end

  sort!(allTimes)
  allTimes = unique(allTimes)

  isPresent = similar(times, BitArray{1})

  for i in eachindex(times)
    isPresent[i] = BitArray{1}(length(times[i]))
    for j in eachindex(allTimes)
      isPresent[i][j] = findfirst(times[i], allTimes[j]) > 0 ? true : false
    end
  end

  return allTimes, isPresent
end

function is_in_subset(mainSet::Vector{Float64}, subSet::Vector{Float64})
  result = falses(length(mainSet))
  dimsubSet = length(subSet)
  if dimsubSet == 0
    return result
  end
  dimSet = length(mainSet)

  setElement = subSetElement = 0.0

  dimSet >= dimsubSet || error("set is required to be larger or equal than subset")

  for i in eachindex(mainSet) # loop in set
    j = 1
    setElement = mainSet[i]
    while true # loop in subset
      subSetElement = subSet[j]
      result[i] = false
      # if smaller no hope, leave false and go to next i
      if setElement < subSetElement
        break
      end
      # match!
      if setElement == subSetElement
        result[i] = true
        break
      end

      # if larger, leave false if at end or go to next j
      if j == dimsubSet
        break
      end

      j += 1
    end
  end

  return result
end
