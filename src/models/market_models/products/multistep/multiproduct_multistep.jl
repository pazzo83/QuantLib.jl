type MultiProductMultiStepCommon
  rateTimes::Vector{Float64}
  evolution::EvolutionDescription
end

function MultiProductMultiStepCommon(rateTimes::Vector{Float64})
  n = length(rateTimes) - 1
  evolutionTimes = zeros(n)
  relevanceRates = Vector{Pair{Int, Int}}(n)
  @simd for i in eachindex(evolutionTimes)
    @inbounds evolutionTimes[i] = rateTimes[i]
    @inbounds relevanceRates[i] = Pair(i, i+1)
  end

  evolution = EvolutionDescription(rateTimes, evolutionTimes, relevanceRates)

  return MultiProductMultiStepCommon(rateTimes, evolution)
end

clone(mpmsc::MultiProductMultiStepCommon) = MultiProductMultiStepCommon(copy(mpmsc.rateTimes), clone(mpmsc.evolution))

rate_times(mpms::MultiProductMultiStep) = mpms.common.rateTimes
get_evolution(mpms::MultiProductMultiStep) = mpms.common.evolution

type MultiStepInverseFloater <: MultiProductMultiStep
  common::MultiProductMultiStepCommon
  fixedAccruals::Vector{Float64}
  floatingAccruals::Vector{Float64}
  fixedStrikes::Vector{Float64}
  fixedMultipliers::Vector{Float64}
  floatingSpreads::Vector{Float64}
  paymentTimes::Vector{Float64}
  payer::Bool
  multiplier::Float64
  lastIndex::Int
  currentIndex::Int
end

function MultiStepInverseFloater(rateTimes::Vector{Float64},
                                fixedAccruals::Vector{Float64},
                                floatingAccruals::Vector{Float64},
                                fixedStrikes::Vector{Float64},
                                fixedMultipliers::Vector{Float64},
                                floatingSpreads::Vector{Float64},
                                paymentTimes::Vector{Float64},
                                payer::Bool = true)

  common = MultiProductMultiStepCommon(rateTimes)
  multiplier = payer ? -1.0 : 1.0

  return MultiStepInverseFloater(common, fixedAccruals, floatingAccruals, fixedStrikes, fixedMultipliers, floatingSpreads, paymentTimes,
                                payer, multiplier, length(rateTimes), 1)
end

type ExerciseAdapter{M <: MarketModelExerciseValue} <: MultiProductMultiStep
  common::MultiProductMultiStepCommon
  exercise::M
  numberOfProducts::Int
  isExerciseTime::BitArray{1}
  currentIndex::Int
end

function ExerciseAdapter(exercise::MarketModelExerciseValue, numberOfProducts::Int = 1)
  ex = clone(exercise)

  common = MultiProductMultiStepCommon(ex.evolution.rateTimes)
  isExerciseTime = ex.isExerciseTime

  return ExerciseAdapter(common, ex, numberOfProducts, isExerciseTime, 1)
end

number_of_products(::MultiStepInverseFloater) = 1

possible_cash_flow_times(msif::MultiStepInverseFloater) = msif.paymentTimes
possible_cash_flow_times(ea::ExerciseAdapter) = possible_cash_flow_times(ea.exercise)

max_number_of_cashflows_per_product_per_step(::ExerciseAdapter) = 1
max_number_of_cashflows_per_product_per_step(::MultiStepInverseFloater) = 1

reset!(msif::MultiStepInverseFloater) = msif.currentIndex = 1

function reset!(ea::ExerciseAdapter)
  reset!(ea.exercise)
  ea.currentIndex = 1

  return ea
end

function next_time_step!(msif::MultiStepInverseFloater, currentState::CurveState, numberCashFlowsThisStep::Vector{Int}, genCashFlows::Vector{Vector{MarketModelCashFlow}})
  liborRate = forward_rate(currentState, msif.currentIndex)
  inverseFloatingCoupon = max((msif.fixedStrikes[msif.currentIndex] - msif.fixedMultipliers[msif.currentIndex] * liborRate), 0.0) * msif.fixedAccruals[msif.currentIndex]
  floatingCoupon = (liborRate + msif.floatingSpreads[msif.currentIndex]) * msif.floatingAccruals[msif.currentIndex]

  genCashFlows[1][1].timeIndex = msif.currentIndex
  genCashFlows[1][1].amount = msif.multiplier * (inverseFloatingCoupon - floatingCoupon)

  numberCashFlowsThisStep[1] = 1
  msif.currentIndex += 1

  return (msif.currentIndex == msif.lastIndex)
end

function next_time_step!(ea::ExerciseAdapter, currentState::CurveState, numberCashFlowsThisStep::Vector{Int}, genCashFlows::Vector{Vector{MarketModelCashFlow}})
  fill!(numberCashFlowsThisStep, 0)
  isDone = false

  next_step!(ea.exercise, currentState)
  if ea.isExerciseTime[ea.currentIndex]
    cf = get_value(ea.exercise, currentState)
    numberCashFlowsThisStep[1] = 1
    genCashFlows[1][1] = cf
    isDone = true
  end

  ea.currentIndex += 1
  return isDone || ea.currentIndex == length(ea.isExerciseTime)+1
end

## Clone ##
clone(msif::MultiStepInverseFloater) = MultiStepInverseFloater(clone(msif.common), copy(msif.fixedAccruals), copy(msif.floatingAccruals), copy(msif.fixedStrikes),
                                        copy(msif.fixedMultipliers), copy(msif.floatingSpreads), copy(msif.paymentTimes), msif.payer, msif.multiplier,
                                        msif.lastIndex, msif.currentIndex)

clone(ea::ExerciseAdapter) = ExerciseAdapter{typeof(ea.exercise)}(clone(ea.common), clone(ea.exercise), ea.numberOfProducts, ea.isExerciseTime, ea.currentIndex)
