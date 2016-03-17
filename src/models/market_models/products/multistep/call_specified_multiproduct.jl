type CallSpecifiedMultiProduct{M <: MarketModelMultiProduct, E <: ExerciseStrategy, M2 <: MarketModelMultiProduct}
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
  n = max_number_of_cashflows_per_step(rebate)

  dummyCashFlowsGenerated = fill(Vector{MarketModelCashFlow}(n), products)

  return CallSpecifiedMultiProduct(underlying, strategy, rebate, evolution, isPresent, cashFlowTimes, rebateOffset, false, dummyCashFlowsThisStep,
                                  dummyCashFlowsGenerated, -1, true)
end
