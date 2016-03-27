type CallSpecifiedPathwiseMultiProduct{M <: MarketModelPathwiseMultiProduct, E <: ExerciseStrategy, M2 <: MarketModelPathwiseMultiProduct} <: MarketModelPathwiseMultiProduct
  underlying::M
  strategy::E
  rebate::M2
  evolution::EvolutionDescription
  isPresent::Vector{BitArray{1}}
  cashFlowTimes::Vector{Float64}
  rebateOffset::Int
  wasCalled::Bool
  dummyCashFlowsThisStep::Vector{Int}
  dummyCashFlowsGenerated::Vector{Vector{MarketModelPathWiseCashFlow}}
  currentIndex::Int
  callable::Bool
end

function CallSpecifiedPathwiseMultiProduct(underlying::MarketModelPathwiseMultiProduct,
                                  strategy::ExerciseStrategy)
  underlying = clone(underlying)
  strategy = clone(strategy)

  products = number_of_products(underlying)
  d1 = get_evolution(underlying)
  rateTimes1 = d1.rateTimes
  evolutionTimes1 = d1.evolutionTimes
  exerciseTimes = strategy.exerciseTimes

  # TODO check if rebate is there
  description = EvolutionDescription(rateTimes1, exerciseTimes)
  amounts = zeros(products, length(exerciseTimes))
  rebate = MarketModelPathwiseCashRebate(description, exerciseTimes, amounts, products)

  allEvolutionTimes = Vector{Vector{Float64}}(4)
  allEvolutionTimes[1] = evolutionTimes1
  allEvolutionTimes[2] = exerciseTimes
  allEvolutionTimes[3] = rebate.evolution.evolutionTimes
  allEvolutionTimes[4] = relevant_times(strategy)

  mergedEvolutionTimes, isPresent = merge_times(allEvolutionTimes)

  # TODO Add relevant rates
  evolution = EvolutionDescription(rateTimes1, mergedEvolutionTimes)
  cashFlowTimes = possible_cash_flow_times(underlying)
  rebateOffset = length(cashFlowTimes)
  rebateTimes = possible_cash_flow_times(rebate)
  cashFlowTimes = vcat(cashFlowTimes, rebateTimes)

  dummyCashFlowsThisStep = zeros(Int, products)
  n = max_number_of_cashflows_per_product_per_step(rebate)

  dummyCashFlowsGenerated = fill(MarketModelPathWiseCashFlow[MarketModelPathWiseCashFlow(d1.numberOfRates + 1) for i = 1:n], products)

  return CallSpecifiedPathwiseMultiProduct(underlying, strategy, rebate, evolution, isPresent, cashFlowTimes, rebateOffset, false, dummyCashFlowsThisStep,
                                  dummyCashFlowsGenerated, 1, true)
end

number_of_products(cs::CallSpecifiedPathwiseMultiProduct) = number_of_products(cs.underlying)
already_deflated(cs::CallSpecifiedPathwiseMultiProduct) = already_deflated(cs.underlying)
possible_cash_flow_times(cs::CallSpecifiedPathwiseMultiProduct) = cs.cashFlowTimes
max_number_of_cashflows_per_product_per_step(cs::CallSpecifiedPathwiseMultiProduct) =
            max(max_number_of_cashflows_per_product_per_step(cs.underlying), max_number_of_cashflows_per_product_per_step(cs.rebate))

clone(cs::CallSpecifiedPathwiseMultiProduct) = CallSpecifiedPathwiseMultiProduct(clone(cs.underlying), clone(cs.strategy), clone(cs.rebate), clone(cs.evolution), deepcopy(cs.isPresent),
                                        copy(cs.cashFlowTimes), cs.rebateOffset, cs.wasCalled, copy(cs.dummyCashFlowsThisStep), deepcopy(cs.dummyCashFlowsGenerated),
                                        cs.currentIndex, cs.callable)
