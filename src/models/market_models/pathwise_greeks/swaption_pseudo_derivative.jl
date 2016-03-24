type SwaptionPseudoDerivative{M <: AbstractMarketModel}
  inputModel::M
  varianceDerivatives::Vector{Matrix{Float64}}
  volatilityDerivatives::Vector{Matrix{Float64}}
  impliedVolatility::Float64
  expiry::Float64
  variance::Float64
end

function SwaptionPseudoDerivative(inputModel::AbstractMarketModel, startIndex::Int, endIndex::Int)
  varianceDerivatives = Vector{Matrix{Float64}}()
  subRateTimes = inputModel.evolution.rateTimes[startIndex:endIndex]

  subForwards = inputModel.initialRates[startIndex:endIndex - 1]

  cs = LMMCurveState(subRateTimes)
  set_on_forward_rates!(cs, subForwards)

  zed = coterminal_swap_zed_matrix(SwapForwardMappings(), cs, inputModel.displacements[1])
  factors = inputModel.numberOfFactors

  # first compute the variance and implied vol
  variance_ = 0.0
  idx = 1

  while idx <= inputModel.evolution.numberOfSteps && inputModel.evolution.firstAliveRate[idx] <= startIndex
    thisPseudo = inputModel.pseudoRoots[idx]

    thisVariance = 0.0
    for j = startIndex:endIndex, k = startIndex:endIndex, f = 1:factors
      thisVariance += zed[1, j - (startIndex - 1)] * thisPseudo[j, f] * thisPseudo[k, f] * zed[1, k - (startIndex - 1)]
    end

    variance_ += thisVariance
    idx += 1
  end

  stopIndex = idx

  expiry = subRateTimes[1]

  impliedVolatility = sqrt(variance_ / expiry)

  scale_ = 0.5 * (1.0 / expiry) / impliedVolatility

  numberRates = inputModel.evolution.numberOfRates

  thisDerivative = zeros(numberRates, factors)
  nullDerivative = zeros(numberRates, factors)

  idx = 1

  while idx <= stopIndex
    thisPseudo = inputModel.pseudoRoots[idx]

    for rate1 = startIndex:endIndex
      zIndex = rate1 - (startIndex - 1)
      for f = 1:factors
        sum_ = 0.0
        for rate2 = startIndex:endIndex
          zIndex2 = rate2 - (startIndex - 1)
          sum_ += zed[1, zIndex2] * thisPseudo[rate2, f]
        end
        sum_ *= 2.0 * zed[1, zIndex]
        thisDerivative[rate1, f] = sum_
      end
    end

    push!(varianceDerivatives, copy(thisDerivative))

    for rate1 = startIndex:endIndex, f = 1:factors
      thisDerivative[rate1, f]
    end

    push!(volatilityDerivatives, copy(thisDerivative))

    idx += 1
  end

  while idx <= inputModel.evolution.numberOfSteps
    push!(varianceDerivatives, copy(nullDerivative))
    push!(volatilityDerivatives, copy(nullDerivative))
  end

  return SwaptionPseudoDerivative(inputModel, varianceDerivatives, volatilityDerivatives, impliedVolatility, expiry, variance_)
end

type CapPseudoDerivative{M <: AbstractMarketModel}
  inputModel::M
  varianceDerivatives::Vector{Matrix{Float64}}
  volatilityDerivatives::Vector{Matrix{Float64}}
  impliedVolatility::Float64
  vega::Float64
  firstDF::Float64
end
