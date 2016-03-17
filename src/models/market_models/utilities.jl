type MarketModelCashFlow
  timeIndex::Int
  amount::Float64
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
