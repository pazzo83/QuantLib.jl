function flat_vol_covariance(t1::Float64, t2::Float64, T::Float64, S::Float64, v1::Float64, v2::Float64)
  t1 <= t2 || error("integrations bounds are in reverse order")
  cutOff = min(S, T)
  if t1 >= cutOff
    return 0.0
  else
    cutOff = min(t2, cutOff)
    return (cutOff - t1) * v1 * v2
  end
end

struct FlatVol <: AbstractMarketModel
  covariance::Vector{Matrix{Float64}}
  totalCovariance::Vector{Matrix{Float64}}
  numberOfFactors::Int
  numberOfRates::Int
  numberOfSteps::Int
  initialRates::Vector{Float64}
  displacements::Vector{Float64}
  evolution::EvolutionDescription
  pseudoRoots::Vector{Matrix{Float64}}
end

function FlatVol(vols::Vector{Float64},
                corr::PiecewiseConstantCorrelation,
                evolution::EvolutionDescription,
                numberOfFactors::Int,
                initialRates::Vector{Float64},
                displacements::Vector{Float64})

  numberOfRates = length(initialRates)
  numberOfSteps = length(evolution.evolutionTimes)
  pseudoRoots = Matrix{Float64}[Matrix{Float64}(undef, numberOfRates, numberOfFactors) for i = 1:numberOfSteps]

  rateTimes = evolution.rateTimes
  numberOfRates == length(rateTimes) - 1 || error("mismatch between number of rates and rate times")
  numberOfRates == length(displacements) || error("mismatch between number of rates and displacements")
  numberOfRates == length(vols) || error("mismatch between number of rates and vols")
  numberOfRates <= numberOfFactors * numberOfSteps || error("number of rates greater than number of factors times number of steps")
  numberOfFactors <= numberOfRates || error("number of factors cannot be greater than number of rates")
  numberOfFactors > 0 || error("number of factors must be greater than zero")

  effStopTime = 0.0
  corrTimes = corr.times
  evolTimes = evolution.evolutionTimes
  covariance = Matrix{Float64}(undef, numberOfRates, numberOfRates)

  kk = 1
  @inbounds @simd for k = 1:numberOfSteps
    # one covariance per evolution step
    fill!(covariance, 0.0)

    # there might be more than one correlation matrix in a single evolution step
    while corrTimes[kk] < evolTimes[k]
      effStartTime = effStopTime
      corrMatrix = corr.correlations[kk]
      for i = 1:numberOfRates, j = i:numberOfRates
        cov = flat_vol_covariance(effStartTime, effStopTime, rateTimes[i], rateTimes[j], vols[i], vols[j])
        covariance[i, j] += cov * corrMatrix[i, j]
      end

      kk += 1
    end

    # last part in the evolution step
    effStartTime = effStopTime
    effStopTime = evolTimes[k]
    corrMatrix = corr.correlations[kk]
    for i = 1:numberOfRates, j = 1:numberOfRates
      cov = flat_vol_covariance(effStartTime, effStopTime, rateTimes[i], rateTimes[j], vols[i], vols[j])
      covariance[i, j] += cov * corrMatrix[i, j]
    end

    # no more use for the kk-th correlation matrix
    while kk < length(corrTimes) && corrTimes[kk] <= evolTimes[k]
      kk+=1
    end

    # make it symmetric
    for i = 1:numberOfRates, j = i+1:numberOfRates
      covariance[j, i] = covariance[i, j]
    end
    pseudoRoots[k] = rank_reduced_sqrt(covariance, numberOfFactors, 1.0, NoneSalvagingAlgo())
  end
  return FlatVol(Vector{Matrix{Float64}}(), Vector{Matrix{Float64}}(), numberOfFactors, numberOfRates, numberOfSteps, initialRates, displacements, evolution, pseudoRoots)
end

clone(flatVol::FlatVol) = FlatVol(deepcopy(flatVol.covariance), deepcopy(flatVol.totalCovariance), flatVol.numberOfFactors, flatVol.numberOfRates, flatVol.numberOfSteps,
                                  copy(flatVol.initialRates), copy(flatVol.displacements), clone(flatVol.evolution), deepcopy(flatVol.pseudoRoots))
