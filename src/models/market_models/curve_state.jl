type LMMCurveState <: CurveState
  numberOfRates::Int
  rateTimes::Vector{Float64}
  rateTaus::Vector{Float64}
  firstIdx::Int
  discRatios::Vector{Float64}
  forwardRates::Vector{Float64}
  cmSwapRates::Vector{Float64}
  cmSwapAnnuties::Vector{Float64}
  cotSwapRates::Vector{Float64}
  cotAnnuities::Vector{Float64}
  firstCotAnnuityComped::Int
end

function LMMCurveState(rateTimes::Vector{Float64})
  numberOfRates = isempty(rateTimes) ? 0 : length(rateTimes) - 1
  rateTaus = check_increasing_times_and_calculate_taus(rateTimes)

  discRatios = ones(numberOfRates + 1)
  forwardRates = Vector{Float64}(numberOfRates)
  cmSwapRates = Vector{Float64}(numberOfRates)
  cmSwapAnnuties = fill(rateTaus[numberOfRates], numberOfRates)
  cotSwapRates = Vector{Float64}(numberOfRates)
  cotAnnuities = fill(rateTaus[numberOfRates], numberOfRates)

  return LMMCurveState(numberOfRates, rateTimes, rateTaus, numberOfRates+1, discRatios, forwardRates, cmSwapRates, cmSwapAnnuties, cotSwapRates,
                      cotAnnuities, numberOfRates+1)
end

function set_on_forward_rates!(lmm::LMMCurveState, rates::Vector{Float64}, firstValidIndex::Int=1)
  length(rates) == lmm.numberOfRates || error("rates mismatch")
  firstValidIndex <= lmm.numberOfRates || error("first valid index must be less than or equal to numberOfRates")

  # first copy input
  lmm.firstIdx = firstValidIndex
  lmm.forwardRates[lmm.firstIdx:end] = rates[lmm.firstIdx:end]

  # then calculate the discount ratios
  # taken care at constructor time discRatios[numberOfRates] = 1.0
  for i = lmm.firstIdx:lmm.numberOfRates
    lmm.discRatios[i+1] = lmm.discRatios[i] / (1.0 + lmm.forwardRates[i] * lmm.rateTaus[i])
  end

  # lazy evaluation of :
  # - coterminal swap rates/annuities
  # - constant maturity swap rates/annuities

  lmm.firstCotAnnuityComped = lmm.numberOfRates+1
end

function discount_ratio(lmm::LMMCurveState, i::Int, j::Int)
  lmm.firstIdx <= lmm.numberOfRates || error("curve state not initialized yet")
  min(i, j) >= lmm.firstIdx || error("invalid index")
  max(i, j) <= lmm.numberOfRates + 1 || error("invalid index")

  return lmm.discRatios[i] / lmm.discRatios[j]
end

function coterminal_from_discount_ratios!(lmm::LMMCurveState, firstValidIndex::Int, discountFactors::Vector{Float64}, taus::Vector{Float64},
                                        cotSwapRates::Vector{Float64}, cotSwapAnnuities::Vector{Float64})
  nCotSwapRates = length(cotSwapRates)
  length(taus) == nCotSwapRates || error("taus size != cotSwapRate size")
  length(cotSwapAnnuities) == nCotSwapRates || error("cotSwapAnnuities size != cotSwapRate size")
  length(discountFactors) == nCotSwapRates + 1 || error("discountFactors size != cotSwapRate size + 1")

  cotSwapAnnuities[nCotSwapRates] = taus[nCotSwapRates] * discountFactors[nCotSwapRates + 1]
  cotSwapRates[nCotSwapRates] = (discountFactors[nCotSwapRates] - discountFactors[nCotSwapRates + 1]) / cotSwapAnnuities[nCotSwapRates]

  for i = nCotSwapRates:-1:firstValidIndex + 1
    cotSwapAnnuities[i-1] = cotSwapAnnuities[i] + taus[i-1] * discountFactors[i]
    cotSwapRates[i-1] = (discountFactors[i-1] - discountFactors[nCotSwapRates+1]) / cotSwapAnnuities[i-1]
  end

  return cotSwapAnnuities, cotSwapRates
end

function forward_rate(lmm::LMMCurveState, i::Int)
  lmm.firstIdx <= lmm.numberOfRates || error("curve state not initialized yet")
  (i >= lmm.firstIdx && i <= lmm.numberOfRates) || error("invalid index")

  return lmm.forwardRates[i]
end

function coterminal_swap_annuity(lmm::LMMCurveState, numeraire::Int, i::Int)
  lmm.firstIdx <= lmm.numberOfRates || error("curve state not initialized yet")
  (numeraire >= lmm.firstIdx && numeraire <= lmm.numberOfRates + 1) || error("invalid numeraire")
  (i >= lmm.firstIdx && i <= lmm.numberOfRates) || error("invalid index")

  if lmm.firstCotAnnuityComped <= i
    return lmm.cotAnnuities[i] / lmm.discRatios[numeraire]
  end

  if lmm.firstCotAnnuityComped == lmm.numberOfRates+1
    lmm.cotAnnuities[lmm.numberOfRates] = lmm.rateTaus[lmm.numberOfRates] * lmm.discRatios[lmm.numberOfRates + 1]
    lmm.firstCotAnnuityComped -= 1
  end

  for j = lmm.firstCotAnnuityComped-1:-1:i
    lmm.cotAnnuities[j] = lmm.cotAnnuities[j+1] + lmm.rateTaus[j] * lmm.discRatios[j+1]
  end

  lmm.firstCotAnnuityComped = i

  return lmm.cotAnnuities[i] / lmm.discRatios[numeraire]
end

function coterminal_swap_rate(lmm::LMMCurveState, i::Int)
  lmm.firstIdx <= lmm.numberOfRates || error("curve state not initialized yet")
  (i >= lmm.firstIdx && i <= lmm.numberOfRates) || error("invalid index")

  res = (lmm.discRatios[i] / lmm.discRatios[lmm.numberOfRates+1] - 1.0) / coterminal_swap_annuity(lmm, lmm.numberOfRates+1, i)

  return res
end

function coterminal_swap_rates!(lmm::LMMCurveState)
  lmm.firstIdx <= lmm.numberOfRates || error("curve state not initialized yet")
  coterminal_from_discount_ratios!(lmm, lmm.firstIdx, lmm.discRatios, lmm.rateTaus, lmm.cotSwapRates, lmm.cotAnnuities)

  return lmm.cotSwapRates
end

clone(lmm::LMMCurveState) = LMMCurveState(lmm.numberOfRates, copy(lmm.rateTimes), copy(lmm.rateTaus), lmm.firstIdx, copy(lmm.discRatios), copy(lmm.forwardRates),
                                          copy(lmm.cmSwapRates), copy(lmm.cmSwapAnnuties), copy(lmm.cotSwapRates), copy(lmm.cotAnnuities), lmm.firstCotAnnuityComped)
