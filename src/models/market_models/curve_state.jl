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

  return LMMCurveState(numberOfRates, rateTimes, rateTaus, numberOfRates, discRatios, forwardRates, cmSwapRates, cmSwapAnnuties, cotSwapRates,
                      cotAnnuities, numberOfRates)
end
