type CallSpecifiedMultiProduct{M <: MarketModelMultiProduct, E <: ExerciseStrategy, M2 <: MarketModelMultiProduct} <: MarketModelMultiProduct
  underlying::M
  strategy::E
  rebate::M2
  evolution::EvolutionDescription
  isPresent::Vector{BitArray{1}}
  cashFlowTimes::Vector{Float64}
  rebateOffset::Int
  wasCalled::Bool
  dummyCashFlowsThisStep::Vector{Int}
  dummyCashFlowsGenerated::Vector{Vector{MarketModelCashFlow}}
  currentIndex::Int
  callable::Bool
end

function CallSpecifiedMultiProduct(underlying::MarketModelMultiProduct,
                                  strategy::ExerciseStrategy,
                                  rebate::MarketModelMultiProduct)

  underlying = clone(underlying)
  strategy = clone(strategy)
  rebate = clone(rebate)

  products = number_of_products(underlying)
  d1 = get_evolution(underlying)
  rateTimes1 = d1.rateTimes
  evolutionTimes1 = d1.evolutionTimes
  exerciseTimes = strategy.exerciseTimes

  # TODO check if rebate isn't there
  d2 = get_evolution(rebate)
  rateTimes2 = d2.rateTimes
  rateTimes1 == rateTimes2 || error("incompatible rate times")

  allEvolutionTimes = Vector{Vector{Float64}}(4)
  allEvolutionTimes[1] = evolutionTimes1
  allEvolutionTimes[2] = exerciseTimes
  allEvolutionTimes[3] = d2.evolutionTimes
  allEvolutionTimes[4] = relevant_times(strategy)

  mergedEvolutionTimes, isPresent = merge_times(allEvolutionTimes)

  evolution = EvolutionDescription(rateTimes1, mergedEvolutionTimes)
  cashFlowTimes = possible_cash_flow_times(underlying)
  rebateOffset = length(cashFlowTimes)
  rebateTimes = possible_cash_flow_times(rebate)
  cashFlowTimes = vcat(cashFlowTimes, rebateTimes)

  dummyCashFlowsThisStep = zeros(Int, products)
  n = max_number_of_cashflows_per_product_per_step(rebate)

  dummyCashFlowsGenerated = fill(Vector{MarketModelCashFlow}(n), products)

  return CallSpecifiedMultiProduct(underlying, strategy, rebate, evolution, isPresent, cashFlowTimes, rebateOffset, false, dummyCashFlowsThisStep,
                                  dummyCashFlowsGenerated, 1, true)
end

get_evolution(cs::CallSpecifiedMultiProduct) = cs.evolution
number_of_products(cs::CallSpecifiedMultiProduct) = number_of_products(cs.underlying)
possible_cash_flow_times(cs::CallSpecifiedMultiProduct) = cs.cashFlowTimes
function reset!(cs::CallSpecifiedMultiProduct)
  reset!(cs.underlying)
  reset!(cs.rebate)
  reset!(cs.strategy)
  cs.currentIndex = 1
  cs.wasCalled = false

  return cs
end

max_number_of_cashflows_per_product_per_step(cs::CallSpecifiedMultiProduct) =
              max(max_number_of_cashflows_per_product_per_step(cs.underlying), max_number_of_cashflows_per_product_per_step(cs.rebate))

function next_time_step!(cs::CallSpecifiedMultiProduct, currentState::CurveState, numberCashFlowsThisStep::Vector{Int}, genCashFlows::Vector{Vector{MarketModelCashFlow}})
  isUnderlyingTime = cs.isPresent[1][cs.currentIndex]
  isExerciseTime = cs.isPresent[2][cs.currentIndex]
  isRebateTime = cs.isPresent[3][cs.currentIndex]
  isStrategyRelevantTime = cs.isPresent[4][cs.currentIndex]

  isDone = false

  if ~cs.wasCalled && isStrategyRelevantTime
    next_step!(cs.strategy, currentState)
  end

  if ~cs.wasCalled && isExerciseTime && cs.callable
    cs.wasCalled = get_exercise!(cs.strategy, currentState)
  end

  if cs.wasCalled
    if isRebateTime
      isDone = next_time_step!(cs.rebate, currentState, numberCashFlowsThisStep, genCashFlows)
      @inbounds @simd for i in eachindex(numberCashFlowsThisStep)
        for j = 1:numberCashFlowsThisStep[i]
          genCashFlows[i][j].timeIndex += cs.rebateOffset
        end
      end
    end
  else
    if isRebateTime
      next_time_step!(cs.rebate, currentState, cs.dummyCashFlowsThisStep, cs.dummyCashFlowsGenerated)
    end

    if isUnderlyingTime
      isDone = next_time_step!(cs.underlying, currentState, numberCashFlowsThisStep, genCashFlows)
    end
  end

  cs.currentIndex += 1
  return isDone || cs.currentIndex == length(cs.evolution.evolutionTimes) + 1
end

clone(cs::CallSpecifiedMultiProduct) = CallSpecifiedMultiProduct(clone(cs.underlying), clone(cs.strategy), clone(cs.rebate), clone(cs.evolution), deepcopy(cs.isPresent),
                                        copy(cs.cashFlowTimes), cs.rebateOffset, cs.wasCalled, copy(cs.dummyCashFlowsThisStep), deepcopy(cs.dummyCashFlowsGenerated),
                                        cs.currentIndex, cs.callable)
