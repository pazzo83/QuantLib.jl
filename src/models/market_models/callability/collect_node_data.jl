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
  #cashFlowsGenerated[1] = Vector{MarketModelCashFlow}(max_number_of_cashflows_per_product_per_step(product))
  cashFlowsGenerated[1] = MarketModelCashFlow[MarketModelCashFlow() for i = 1:max_number_of_cashflows_per_product_per_step(product)]

  rateTimes = get_evolution(product).rateTimes

  cashFlowTimes = possible_cash_flow_times(product)
  rebateTimes = possible_cash_flow_times(rebate)
  controlTimes = possible_cash_flow_times(control)

  n = length(cashFlowTimes)
  productDiscounters = Vector{MarketModelDiscounter}(n)
  @simd for i in eachindex(productDiscounters)
    @inbounds productDiscounters[i] = MarketModelDiscounter(cashFlowTimes[i], rateTimes)
  end

  n = length(rebateTimes)
  rebateDiscounters = Vector{MarketModelDiscounter}(n)
  @simd for i in eachindex(rebateDiscounters)
    @inbounds rebateDiscounters[i] = MarketModelDiscounter(rebateTimes[i], rateTimes)
  end

  n = length(controlTimes)
  controlDiscounters = Vector{MarketModelDiscounter}(n)
  @simd for i in eachindex(controlDiscounters)
    @inbounds controlDiscounters[i] = MarketModelDiscounter(controlTimes[i], rateTimes)
  end

  evolution = get_evolution(product)
  numeraires = evolver.numeraires

  evolutionTimes = evolution.evolutionTimes

  isProductTime = is_in_subset(evolutionTimes, get_evolution(product).evolutionTimes)
  isRebateTime = is_in_subset(evolutionTimes, rebate.evolution.evolutionTimes)
  isControlTime = is_in_subset(evolutionTimes, control.evolution.evolutionTimes)
  isBasisTime = is_in_subset(evolutionTimes, dataProvider.evolution.evolutionTimes)
  isExerciseTime = falses(length(evolutionTimes))
  v = rebate.isExerciseTime
  exercises = idx = 1

  @inbounds @simd for i in eachindex(evolutionTimes)
    if isRebateTime[i]
      if v[idx]
        isExerciseTime[i] = true
        exercises += 1
      end
      idx += 1
    end
  end

  # collectedData = Vector{Vector{NodeData}}(exercises+1)
  resize!(collectedData, exercises)
  for i in eachindex(collectedData)
    # collectedData[i] = Vector{NodeData}(numberOfPaths)
    collectedData[i] = NodeData[NodeData() for i = 1:numberOfPaths]
  end

  for i = 1:numberOfPaths
    start_new_path!(evolver)
    reset!(product)
    reset!(rebate)
    reset!(control)
    reset!(dataProvider)

    principalInNumerairePortfolio = 1.0
    isDone = false
    nextExercise = 1
    collectedData[1][i].cumulatedCashFlows = 0.0

    while true
      currentStep = evolver.currentStep
      advance_step!(evolver)
      currentState = current_state(evolver)
      numeraire = numeraires[currentStep]

      if isRebateTime[currentStep]
        next_step!(rebate, currentState)
      end

      if isControlTime[currentStep]
        next_step!(control, currentState)
      end

      if isBasisTime[currentStep]
        next_step!(dataProvider, currentState)
      end

      if isExerciseTime[currentStep]
        data = collectedData[nextExercise+1][i]

        exerciseValue = get_value(rebate, currentState)
        data.exerciseValue = exerciseValue.amount * numeraire_bonds(rebateDiscounters[exerciseValue.timeIndex], currentState, numeraire) /
                            principalInNumerairePortfolio
        set_values!(dataProvider, currentState, data.values)

        controlValue = get_value(control, currentState)
        data.controlValue = controlValue.amount * numeraire_bonds(controlDiscounters[controlValue.timeIndex], currentState, numeraire) /
                            principalInNumerairePortfolio

        data.cumulatedCashFlows = 0.0
        data.isValid = true
        nextExercise += 1
      end
      if isProductTime[currentStep]
        isDone = next_time_step!(product, currentState, numberCashFlowsThisStep, cashFlowsGenerated)

        for j = 1:numberCashFlowsThisStep[1]
          cf = cashFlowsGenerated[1][j]
          collectedData[nextExercise][i].cumulatedCashFlows += cf.amount * numeraire_bonds(productDiscounters[cf.timeIndex], currentState, numeraire) / principalInNumerairePortfolio
        end
      end

      if ~isDone
        nextNumeraire = numeraires[currentStep + 1]
        principalInNumerairePortfolio *= discount_ratio(currentState, numeraire, nextNumeraire)
      else
        break
      end
    end

    # fill the remaining (un)collected data with nulls
    for j = nextExercise:exercises-1
      data = collectedData[j+1][i]
      data.exerciseValue = data.controlValue = 0.0
      data.cumulatedCashFlows = 0.0
      data.isValid = false
    end
  end
  evolver, collectedData
end
