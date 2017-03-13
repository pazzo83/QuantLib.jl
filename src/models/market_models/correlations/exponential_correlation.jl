function exponential_correlations(rateTimes::Vector{Float64}, longTermCorr::Float64, beta::Float64, gamma::Float64, time_::Float64)
  check_increasing_times(rateTimes)
  (longTermCorr <= 1.0 && longTermCorr >= 0.0) || error("ong term correlation outside of [0,1] interval")
  beta >= 0.0 || error("beta must be greater than or equal to zero")
  (gamma <= 1.0 && gamma >= 0.0) || error("gamma outside [0,1] interval")

  # Calcualte correlation matrix
  nbRows = length(rateTimes)-1
  correlations = zeros(nbRows, nbRows)

  @inbounds @simd for i = 1:nbRows
    # Correlation is defined only between (alive) stochastic rates
    if time_ <= rateTimes[i]
      correlations[i, i] = 1.0
      for j = 1:i
        if time_ <= rateTimes[j]
          correlations[i,j] = correlations[j, i] = longTermCorr + (1.0 - longTermCorr) *
                              exp(-beta * abs(^(rateTimes[i] - time_, gamma) - ^(rateTimes[j] - time_, gamma)))
        end
      end
    end
  end

  return correlations
end

type ExponentialForwardCorrelation <: PiecewiseConstantCorrelation
  numberOfRates::Int
  longTermCorr::Float64
  beta::Float64
  gamma::Float64
  rateTimes::Vector{Float64}
  times::Vector{Float64}
  correlations::Vector{Matrix{Float64}}
end

function ExponentialForwardCorrelation(rateTimes::Vector{Float64},
                                      longTermCorr::Float64 = 0.5,
                                      beta::Float64 = 0.2,
                                      gamma::Float64 = 1.0,
                                      times::Vector{Float64} = Vector{Float64}())
  # build
  numberOfRates = isempty(rateTimes) ? 0 : length(rateTimes) - 1
  numberOfRates > 1 || error("rate times must contain at least two values")

  check_increasing_times(rateTimes)

  # corrTimes must include all rateTimes but the last
  if isempty(times)
    times = rateTimes[1:end-1]
  else
    check_increasing_times(times)
  end

  if is_close(gamma, 1.0)
    temp = rateTimes[1:end-1]
    times == temp || error("corr times must be equal to all rate times (but the last)")

    c = exponential_correlations(rateTimes, longTermCorr, beta, 1.0, 0.0)
    correlations = evolved_matrices(TimeHomogeneousForwardCorrelation, c)
  else
    times[end] <= rateTimes[numberOfRates+1] || error("last corr time is after next-to_last rate time")

    correlations = similar(times, Matrix{Float64})
    time_ = times[1] / 2.0
    correlations[1] = exponential_correlations(rateTimes, longTermCorr, beta, gamma, time_)
    @inbounds @simd for k = 2:length(times)
      time_ = (times[k] + times[k-1]) / 2.0
      correlations[k] = exponential_correlations(rateTimes, longTermCorr, beta, gamma, time)
    end
  end

  return ExponentialForwardCorrelation(numberOfRates, longTermCorr, beta, gamma, rateTimes, times, correlations)
end

clone(efc::ExponentialForwardCorrelation) = ExponentialForwardCorrelation(efc.numberOfRates, efc.longTermCorr, efc.beta, efc.gamma, copy(efc.rateTimes), copy(efc.times), deepcopy(efc.correlations))
