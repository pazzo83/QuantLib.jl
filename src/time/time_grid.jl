type TimeGrid
  times::Vector{Float64}
  dt::Vector{Float64}
  mandatoryTimes::Vector{Float64}
end

function TimeGrid{I <: Integer}(times::Vector{Float64}, steps::I)
  sortedUniqueTimes = sort(unique(times))

  lastTime = sortedUniqueTimes[end]

  # TODO check if steps is 0
  dtMax = lastTime / steps
  periodBegin = 0.0
  times = zeros(1)

  for t in sortedUniqueTimes
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
