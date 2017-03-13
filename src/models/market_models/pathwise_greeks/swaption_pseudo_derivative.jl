type QuickCap
  strike::Float64
  annuities::Vector{Float64}
  currentRates::Vector{Float64}
  expires::Vector{Float64}
  price::Float64
end

function operator(qc::QuickCap)
  function _inner(vol::Float64)
    price = 0.0
    for i in eachindex(qc.annuities)
      price += black_formula(Call(), qc.strike, qc.currentRates[i], vol * sqrt(qc.expires[i]), qc.annuities[i])
    end

    return price - qc.price
  end

  return _inner
end

function get_vega(qc::QuickCap, vol::Float64)
  vega_ = 0.0

  @simd for i in eachindex(qc.annuities)
    @inbounds vega_ += black_formula_vol_derivative(qc.strike, qc.currentRates[i], vol * sqrt(qc.expires[i]), qc.expires[i], qc.annuities[i], 0.0)
  end

  return vega_
end

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
  volatilityDerivatives = Vector{Matrix{Float64}}()
  subRateTimes = inputModel.evolution.rateTimes[startIndex:endIndex+1]

  subForwards = inputModel.initialRates[startIndex:endIndex]

  cs = LMMCurveState(subRateTimes)
  set_on_forward_rates!(cs, subForwards)

  zed = coterminal_swap_zed_matrix(SwapForwardMappings(), cs, inputModel.displacements[1])
  factors = inputModel.numberOfFactors

  # first compute the variance and implied vol
  variance_ = 0.0
  idx = 1

  while idx <= number_of_steps(inputModel.evolution) && inputModel.evolution.firstAliveRate[idx] <= startIndex
    thisPseudo = inputModel.pseudoRoots[idx]

    thisVariance = 0.0
    for j = startIndex:endIndex, k = startIndex:endIndex, f = 1:factors
      @inbounds thisVariance += zed[1, j - (startIndex - 1)] * thisPseudo[j, f] * thisPseudo[k, f] * zed[1, k - (startIndex - 1)]
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

  while idx < stopIndex
    thisPseudo = inputModel.pseudoRoots[idx]

    @inbounds @simd for rate1 = startIndex:endIndex
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
      @inbounds thisDerivative[rate1, f] *= scale_
    end

    push!(volatilityDerivatives, copy(thisDerivative))

    idx += 1
  end

  while idx <= number_of_steps(inputModel.evolution)
    push!(varianceDerivatives, copy(nullDerivative))
    push!(volatilityDerivatives, copy(nullDerivative))
    idx += 1
  end

  return SwaptionPseudoDerivative(inputModel, varianceDerivatives, volatilityDerivatives, impliedVolatility, expiry, variance_)
end

type CapPseudoDerivative{M <: AbstractMarketModel}
  inputModel::M
  volatilityDerivatives::Vector{Matrix{Float64}}
  priceDerivatives::Vector{Matrix{Float64}}
  impliedVolatility::Float64
  vega::Float64
  firstDF::Float64
end

function CapPseudoDerivative(inputModel::AbstractMarketModel,
                            strike::Float64,
                            startIndex::Int,
                            endIndex::Int,
                            firstDF::Float64)

  startIndex <= endIndex || error("for a cap pseudo derivative the start of the cap must be before or equal to the end")
  endIndex <= inputModel.numberOfRates || error("for a cap pseudo derivative the end of the cap must be before the end of the rates")

  priceDerivatives = Vector{Matrix{Float64}}()
  volatilityDerivatives = Vector{Matrix{Float64}}()

  numberCaplets = endIndex - startIndex + 1
  numberRates = inputModel.numberOfRates
  factors = inputModel.numberOfFactors
  curve = LMMCurveState(inputModel.evolution.rateTimes)
  set_on_forward_rates!(curve, inputModel.initialRates)

  totalCovariance = total_covariance(inputModel, inputModel.numberOfSteps)

  displacedImpliedVols = Vector{Float64}(numberCaplets)
  annuities = Vector{Float64}(numberCaplets)
  capletPrices = Vector{Float64}(numberCaplets)
  initialRates = Vector{Float64}(numberCaplets)
  expires = Vector{Float64}(numberCaplets)

  capPrice = 0.0
  guess = 0.0
  minVol = 1e10
  maxVol = 0.0

  @inbounds @simd for j = startIndex:endIndex
    capletIndex = ((j - 1) - (startIndex - 1)) + 1
    resetTime = inputModel.evolution.rateTimes[j]
    expires[capletIndex] = resetTime

    sd = sqrt(totalCovariance[j, j])
    displacedImpliedVols[capletIndex] = sqrt(totalCovariance[j, j] / resetTime)

    forward = inputModel.initialRates[j]
    initialRates[capletIndex] = forward

    annuity = discount_ratio(curve, j+1, 1) * inputModel.evolution.rateTaus[j] * firstDF
    annuities[capletIndex] = forward

    displacement = inputModel.displacements[j]

    guess += displacedImpliedVols[capletIndex] * (forward + displacement) / forward
    minVol = min(minVol, displacedImpliedVols[capletIndex])
    maxVol = max(maxVol, displacedImpliedVols[capletIndex] * (forward + displacement) / forward)

    capletPrice = black_formula(Call(), strike, forward, sd, annuity, displacement)
    capletPrices[capletIndex] = capletPrice
    capPrice += capletPrice
  end

  guess /= numberCaplets

  @inbounds @simd for s = 1:number_of_steps(inputModel.evolution)
    thisDerivative = zeros(numberRates, factors)

    for rate = max(inputModel.evolution.firstAliveRate[s], startIndex):endIndex, f = 1:factors
      expiry = inputModel.evolution.rateTimes[rate]
      volDerivative = inputModel.pseudoRoots[s][rate, f] / (displacedImpliedVols[((rate - 1) - (startIndex - 1) + 1)] * expiry)
      capletVega = black_formula_vol_derivative(strike, inputModel.initialRates[rate], displacedImpliedVols[((rate - 1) - (startIndex - 1) + 1)] * sqrt(expiry), expiry,
                                                annuities[((rate - 1) - (startIndex - 1) + 1)], inputModel.displacements[rate])
      thisDerivative[rate, f] = volDerivative * capletVega
    end
    push!(priceDerivatives, thisDerivative)
  end

  capPricer = QuickCap(strike, annuities, initialRates, expires, capPrice)
  maxEvals = 1000
  accuracy = 1e-6

  solver = BrentSolver(maxEvals)
  impliedVolatility = solve(solver, operator(capPricer), accuracy, guess, minVol * 0.99, maxVol * 1.01)

  vega_ = get_vega(capPricer, impliedVolatility)

  @inbounds @simd for s = 1:number_of_steps(inputModel.evolution)
    thisDerivative = zeros(numberRates, factors)
    for rate = max(inputModel.evolution.firstAliveRate[s], startIndex):endIndex, f = 1:factors
      thisDerivative[rate, f] = priceDerivatives[s][s, f] / vega_
    end

    push!(volatilityDerivatives, thisDerivative)
  end
  return CapPseudoDerivative(inputModel, volatilityDerivatives, priceDerivatives, impliedVolatility, vega_, firstDF)
end
