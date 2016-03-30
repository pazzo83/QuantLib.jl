## General Market Model methods ##
function get_covariance(mm::AbstractMarketModel, i::Int)
  if isempty(mm.covariance)
    resize!(mm.covariance, mm.numberOfSteps)
    for j = 1:mm.numberOfSteps
      mm.covariance[j] = mm.pseudoRoots[j] * mm.pseudoRoots[j]'
    end
    i <= length(mm.covariance) || error("i must be less than or equal to covariance size")
  end

  return mm.covariance[i]
end

function total_covariance(mm::AbstractMarketModel, endIndex::Int)
  if isempty(mm.totalCovariance)
    resize!(mm.totalCovariance, mm.numberOfSteps)
    mm.totalCovariance[1] = get_covariance(mm, 1) # need to trigger calculation
    for j = 2:mm.numberOfSteps
      mm.totalCovariance[j] = mm.totalCovariance[j-1] + mm.covariance[j]
    end
    endIndex <= length(mm.covariance) || error("endIndex must be less than or equal to total covariance size")
  end

  return mm.totalCovariance[endIndex]
end

type MarketModelPathwiseInverseFloater <: MarketModelPathwiseMultiProduct
  rateTimes::Vector{Float64}
  fixedAccruals::Vector{Float64}
  floatingAccruals::Vector{Float64}
  fixedStrikes::Vector{Float64}
  fixedMultipliers::Vector{Float64}
  floatingSpreads::Vector{Float64}
  paymentTimes::Vector{Float64}
  payer::Bool
  multiplier::Float64
  lastIndex::Int
  evolution::EvolutionDescription
  currentIndex::Int
end

function MarketModelPathwiseInverseFloater(rateTimes::Vector{Float64},
                                          fixedAccruals::Vector{Float64},
                                          floatingAccruals::Vector{Float64},
                                          fixedStrikes::Vector{Float64},
                                          fixedMultipliers::Vector{Float64},
                                          floatingSpreads::Vector{Float64},
                                          paymentTimes::Vector{Float64},
                                          payer::Bool)

  multiplier = payer ? -1.0 : 1.0
  lastIndex = length(rateTimes)

  check_increasing_times(paymentTimes)

  length(fixedAccruals) == lastIndex - 1 || error("Incorrect number of fixedAccruals given")
  length(floatingAccruals) == lastIndex - 1 || error("Incorrect number of floatingAccruals given")
  length(fixedStrikes) == lastIndex - 1 || error("Incorrect number of fixedStrikes given")
  length(fixedMultipliers) == lastIndex - 1 || error("Incorrect number of fixedMultipliers given")
  length(floatingSpreads) == lastIndex - 1 || error("Incorrect number of floatingSpreads given")
  length(paymentTimes) == lastIndex - 1 || error("Incorrect number of paymentTimes given")

  evolTimes = copy(rateTimes)
  pop!(evolTimes)

  evolution = EvolutionDescription(rateTimes, evolTimes)

  return MarketModelPathwiseInverseFloater(rateTimes, fixedAccruals, floatingAccruals, fixedStrikes, fixedMultipliers, floatingSpreads,
                                        paymentTimes, payer, multiplier, lastIndex, evolution, 1)
end

number_of_products(::MarketModelPathwiseInverseFloater) = 1
get_evolution(mm::MarketModelPathwiseInverseFloater) = mm.evolution
possible_cash_flow_times(mm::MarketModelPathwiseInverseFloater) = mm.paymentTimes
max_number_of_cashflows_per_product_per_step(mm::MarketModelPathwiseInverseFloater) = 1
already_deflated(::MarketModelPathwiseInverseFloater) = false

reset!(mm::MarketModelPathwiseInverseFloater) = mm.currentIndex = 1

function next_time_step!(mm::MarketModelPathwiseInverseFloater, currentState::CurveState, numberCashFlowsThisStep::Vector{Int},
                        cashFlowsGenerated::Vector{Vector{MarketModelPathWiseCashFlow}})
  numberCashFlowsThisStep[1] = 1
  for i = 2:mm.lastIndex
    cashFlowsGenerated[1][1].amount[i] = 0.0
  end

  liborRate = forward_rate(currentState, mm.currentIndex)
  inverseFloatingCoupon = max((mm.fixedStrikes[mm.currentIndex] - mm.fixedMultipliers[mm.currentIndex] * liborRate), 0.0) * mm.fixedAccruals[mm.currentIndex]
  floatingCoupon = (liborRate + mm.floatingSpreads[mm.currentIndex]) * mm.fixedAccruals[mm.currentIndex]

  cashFlowsGenerated[1][1].timeIndex = mm.currentIndex
  cashFlowsGenerated[1][1].amount[1] = mm.multiplier * (inverseFloatingCoupon - floatingCoupon)

  if inverseFloatingCoupon > 0.0
    cashFlowsGenerated[1][1].amount[mm.currentIndex + 1] = mm.multiplier * ( -mm.fixedMultipliers[mm.currentIndex] * mm.fixedAccruals[mm.currentIndex] -
                                                          mm.floatingAccruals[mm.currentIndex])
  else
    cashFlowsGenerated[1][1].amount[mm.currentIndex + 1] = -mm.multiplier * mm.floatingAccruals[mm.currentIndex]
  end

  mm.currentIndex += 1

  return mm.currentIndex == mm.lastIndex
end

## Clone ##
clone(mm::MarketModelPathwiseInverseFloater) = MarketModelPathwiseInverseFloater(copy(mm.rateTimes), copy(mm.fixedAccruals), copy(mm.floatingAccruals), copy(mm.fixedStrikes),
                                                copy(mm.fixedMultipliers), copy(mm.floatingSpreads), copy(mm.paymentTimes), mm.payer, mm.multiplier,
                                                mm.lastIndex, clone(mm.evolution), mm.currentIndex)
