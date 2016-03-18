function collect_node_data!(evolver::AbstractMarketModelEvolver,
                            product::MarketModelMultiProduct,
                            dataProvider::MarketModelBasisSystem,
                            rebate::MarketModelExerciseValue,
                            control::MarketModelExerciseValue,
                            numberOfPaths::Int,
                            collectedData::Vector{Vector{NodeData}})

  numerairesHeld = Vector{Float64}()
  number_of_products(product) == 1 || error("single product required")

  numberCashFlowsThisStep = Vector{Int}(1)
  cashFlowsGenerated = Vector{Vector{MarketModelCashFlow}}(1)
  cashFlowsGenerated[1] = Vector{MarketModelCashFlow}(max_number_of_cashflows_per_step(product))

  rateTimes = get_evolution(product).rateTimes

  cashFlowTimes = possible_cash_flow_times(product)
  rebateTimes = possible_cash_flow_times(rebate)
  controlTimes = possible_cash_flow_times(control)

  n = length(cashFlowTimes)
  productDiscounters = Vector{MarketModelDiscounter}(n)
  for i in eachindex(productDiscounters)
    productDiscounters[i] = MarketModelDiscounter(cashFlowTimes[i], rateTimes)
  end

  n = length(rebateTimes)
  rebateDiscounters = Vector{MarketModelDiscounter}(n)
  for i in eachindex(rebateDiscounters)
    rebateDiscounters[i] = MarketModelDiscounter(rebateTimes[i], rateTimes))
  end

  n = length(controlTimes)
  controlDiscounters = Vector{MarketModelDiscounter}(n)
  for i in eachindex(controlDiscounters)
    controlDiscounters[i] = MarketModelDiscounter(controlTimes[i], rateTimes)
  end

  evolution = get_evolution(product)
  numeraires = evolver.numeraires

  evolutionTimes = evolution.evolutionTimes

  isProductTime = is_in_subset(evolutionTimes, get_evolution(product).evolutionTimes)
  isRebateTime = is_in_subset(evolutionTimes, rebate.evolution.evolutionTimes)
  isControlTime = is_in_subset(evolutionTimes, control.evolution.evolutionTimes)
  isBasisTimes = is_in_subset(evolutionTimes, dataProvider.evolution.evolutionTimes)
  isExerciseTime = falses(length(evolutionTimes))
  v = rebate.isExerciseTime
  exercises = idx = 1

  for i in eachindex(evolutionTimes)
    if isRebateTime[i]
      if v[idx += 1]
        isExerciseTime[i] = true
        exercise += 1
      end
    end
  end

  collectedData = Vector{Vector{NodeData}}(exercises+1)
  for i in eachindex(collectedData)
    collectedData[i] = Vector{NodeData}(numberOfPaths)
  end
end
