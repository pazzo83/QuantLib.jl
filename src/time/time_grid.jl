using QuantLib

import Base.getindex, Base.lastindex

mutable struct TimeGrid
  times::Vector{Float64}
  dt::Vector{Float64}
  mandatoryTimes::Vector{Float64}
end

function TimeGrid(times::Vector{Float64}, steps::Int)
  sortedUniqueTimes = sort(unique(times))

  lastTime = sortedUniqueTimes[end]

  # TODO check if steps is 0
  dtMax = lastTime / steps
  periodBegin = 0.0
  times = zeros(1)

  @inbounds @simd for t in sortedUniqueTimes
    periodEnd = t
    if periodEnd != 0.0
      nSteps = Int(floor((periodEnd - periodBegin) / dtMax + 0.5))
      nSteps = nSteps != 0 ? nSteps : 1
      dt = (periodEnd - periodBegin) / nSteps

      tempTimes = zeros(nSteps)
      for n=1:nSteps
        tempTimes[n] = periodBegin + n * dt
      end
    end
    periodBegin = periodEnd
    times = vcat(times, tempTimes)
  end

  dt = diff(times)

  return TimeGrid(times, dt, sortedUniqueTimes)
end

function TimeGrid(endTime::Float64, steps::Int)
  endTime > 0.0 || error("negative times not allowed")
  dt = endTime / steps
  times = zeros(steps + 1)
  @simd for i in eachindex(times)
    @inbounds times[i] = dt * (i - 1)
  end

  mandatoryTimes = [endTime]

  dtVec = fill(dt, steps)

  return TimeGrid(times, dtVec, mandatoryTimes)
end

getindex(tg::TimeGrid, i::Int) = tg.times[i]
lastindex(tg::TimeGrid) = lastindex(tg.times)

is_empty(tg::TimeGrid) = length(tg.times) == 0

function closest_index(tg::TimeGrid, t::Float64)
  # stuff
  res = searchsortedfirst(tg.times, t)
  if res == 1
    return 1
  elseif res == length(tg.times) + 1
    return length(tg.times)
  else
    dt1 = tg.times[res] - t
    dt2 = t - tg.times[res - 1]
    if dt1 < dt2
      return res
    else
      return res - 1
    end
  end
end

closest_time(tg::TimeGrid, t::Float64) = tg.times[closest_index(tg, t)]

function return_index(tg::TimeGrid, t::Float64)
  i = closest_index(tg, t)
  if QuantLib.Math.close_enough(t, tg.times[i])
    return i
  else
    error("this time grid is wrong, $i $t $(tg.times[i]) $(tg.times[i-1]) $(tg.times[i+1])")
  end
end
