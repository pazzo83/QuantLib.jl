type MarketModelCashFlow
  timeIndex::Int
  amount::Float64
end

clone(mmcf::MarketModelCashFlow) = MarketModelCashFlow(mmcf.timeIndex, mmcf.amount)

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
