type PathwiseVegasOuterAccountingEngine{P <: MarketModelPathwiseMultiProduct, M <: AbstractMarketModel}
  evolver::LogNormalFwdRateEuler
  product::P
  pseudoRootStructure::M
  vegaBumps::Vector{Vector{Matrix{Float64}}}
  numeraires::Vector{Int}
  initialNumeraireValue::Float64
  numberProducts::Int
  numberRates::Int
  numberCashFlowTimes::Int
  numberSteps::Int
  factors::Int
  numberBumps::Int
  numberElementaryVegas::Int
  jacobianComputers::Vector{RatePseudoRootJacobianAllElements}
  doDeflation::Bool
  currentForwards::Vector{Float64}
  lastForwards::Vector{Float64}
  numerairesHeld::Vector{Float64}
  numberCashFlowsThisStep::Vector{Int}
  cashFlowsGenerated::Vector{Vector{MarketModelPathWiseCashFlow}}
  discounters::Vector{MarketModelPathwiseDiscounter}
  V::Vector{Matrix{Float64}}
  LIBORRatios::Matrix{Float64}
  Discounts::Matrix{Float64}
  StepsDiscountsSquared::Matrix{Float64}
  stepsDiscounts::Vector{Float64}
  LIBORRates::Matrix{Float64}
  partials::Matrix{Float64}
  elementaryVegasThisPath::Vector{Vector{Matrix{Float64}}}
  jacobiansThisPaths::Vector{Vector{Matrix{Float64}}}
  deflatorAndDerivatives::Vector{Float64}
  fullDerivatives::Vector{Float64}
  numberCashFlowsThisIndex::Vector{Vector{Int}}
  totalCashFlowsThisIndex::Vector{Matrix{Float64}}
  cashFlowIndicesThisStep::Vector{Vector{Int}}
end

function PathwiseVegasOuterAccountingEngine(evolver::LogNormalFwdRateEuler,
                                            product::MarketModelPathwiseMultiProduct,
                                            pseudoRootStructure::AbstractMarketModel,
                                            vegaBumps::Vector{Vector{Matrix{Float64}}},
                                            initialNumeraireValue::Float64)
  product = clone(product)
  numberProducts = number_of_products(product)
  doDeflation = ~already_deflated(product)
  numerairesHeld = Vector{Float64}(numberProducts)
  numberCashFlowsThisStep = Vector{Int}(numberProducts)
  cashFlowsGenerated = Vector{Vector{MarketModelPathWiseCashFlow}}(numberProducts)
  stepsDiscounts = Vector{Float64}(pseudoRootStructure.numberOfRates + 1)
  elementaryVegasThisPath = Vector{Vector{Matrix{Float64}}}(numberProducts)
  deflatorAndDerivatives = Vector{Float64}(pseudoRootStructure.numberOfRates + 1)

  stepsDiscounts[1] = 1.0

  numberRates = pseudoRootStructure.numberOfRates
  numberSteps = pseudoRootStructure.numberOfSteps
  factors = pseudoRootStructure.numberOfFactors
  fullDerivatives = Vector{Float64}(numberRates)

  evolution = pseudoRootStructure.evolution
  numeraires = money_market_measure(evolution)

  length(vegaBumps) == numberSteps || error("we need precisely one vector of vega bumps for each step")

  numberBumps = length(vegaBumps[1])

  jacobiansThisPathsModel = Vector{Matrix{Float64}}(numberRates)

  for i in eachindex(jacobiansThisPathsModel)
    jacobiansThisPathsModel[i] = Matrix{Float64}(numberRates, factors)
  end

  jacobianComputers = Vector{RatePseudoRootJacobianAllElements}(numberSteps)
  jacobiansThisPaths = Vector{Vector{Matrix{Float64}}}(numberSteps)
  for i in eachindex(jacobianComputers)
    jacobianComputers[i] = RatePseudoRootJacobianAllElements(pseudoRootStructure.pseudoRoots[i], evolution.firstAliveRate[i],
                            numeraires[i], evolution.rateTaus, pseudoRootStructure.displacements)
    jacobiansThisPaths[i] = deepcopy(jacobiansThisPathsModel)
  end

  VModel = Matrix{Float64}(numberSteps + 1, numberRates)

  Discounts = Matrix{Float64}(numberSteps + 1, numberRates + 1)

  for i = 1:numberSteps+1
    Discounts[i, 1] = 1.0
  end

  V = Vector{Matrix{Float64}}(numberProducts)

  modelCashFlowIndex = Matrix{Float64}(length(possible_cash_flow_times(product)), numberRates + 1)

  numberCashFlowsThisIndex = Vector{Vector{Int}}(numberProducts)
  V = Vector{Matrix{Float64}}(numberProducts)
  totalCashFlowsThisIndex = Vector{Matrix{Float64}}(numberProducts)

  for i in eachindex(numberCashFlowsThisIndex)
    p = max_number_of_cashflows_per_product_per_step(product)
    cashFlowsGenerated[i] = MarketModelPathWiseCashFlow[MarketModelPathWiseCashFlow(numberRates + 1) for i = 1:p]

    numberCashFlowsThisIndex[i] = Vector{Int}(length(possible_cash_flow_times(product)))
    V[i] = copy(VModel)
    totalCashFlowsThisIndex[i] = copy(modelCashFlowIndex)
  end

  LIBORRatios = copy(VModel)
  StepsDiscountsSquared = copy(VModel)
  LIBORRates = copy(VModel)

  cashFlowTimes = possible_cash_flow_times(product)
  numberCashFlowTimes = length(cashFlowTimes)

  rateTimes = product.evolution.rateTimes
  evolutionTimes = product.evolution.evolutionTimes

  discounters = MarketModelPathwiseDiscounter[MarketModelPathwiseDiscounter(cashFlowTimes[j], rateTimes) for j = 1:numberCashFlowTimes]

  ## need to check that we are in money market measure

  # we need to allocate cash-flow times to steps, i.e. what is the last step completed before a flow occurs
  # what we really need for each step, what cash flow time indices to look at

  cashFlowIndicesThisStep = Vector{Int}[Vector{Int}() for i = 1:numberSteps]
  for i in eachindex(cashFlowTimes)
    # idx = upper_bound(evolutionTimes, cashFlowTimes[i])
    # if idx != 1
    #   idx -= 1
    # end
    idx = findlast(evolutionTimes, cashFlowTimes[i])
    if idx == 0
      if cashFlowTimes[i] > evolutionTimes[end]
        idx = length(evolutionTimes)
      else
        idx = 1
      end
    end
    push!(cashFlowIndicesThisStep[idx], i)
  end

  partials = Matrix{Float64}(factors, numberRates)

  begin
    modelVegaMatrix = zeros(numberRates, factors)
    for i in eachindex(elementaryVegasThisPath)
      elementaryVegasThisPath[i] = Matrix[copy(modelVegaMatrix) for j = 1:numberSteps]
    end
  end

  numberElementaryVegas = numberSteps * numberRates * factors

  return PathwiseVegasOuterAccountingEngine(evolver, product, pseudoRootStructure, vegaBumps, numeraires, initialNumeraireValue, numberProducts, numberRates,
                                            numberCashFlowTimes, numberSteps, factors, numberBumps, numberElementaryVegas, jacobianComputers, doDeflation,
                                            Vector{Float64}(), Vector{Float64}(), numerairesHeld, numberCashFlowsThisStep, cashFlowsGenerated, discounters,
                                            V, LIBORRatios, Discounts, StepsDiscountsSquared, stepsDiscounts, LIBORRates, partials, elementaryVegasThisPath,
                                            jacobiansThisPaths, deflatorAndDerivatives, fullDerivatives, numberCashFlowsThisIndex, totalCashFlowsThisIndex,
                                            cashFlowIndicesThisStep)
end

function single_path_values!(pwEng::PathwiseVegasOuterAccountingEngine, values::Vector{Float64})
  initialForwards = copy(pwEng.pseudoRootStructure.initialRates)
  pwEng.currentForwards = copy(initialForwards)

  for i = 1:pwEng.numberProducts
    pwEng.numerairesHeld[i] = 0.0
    for j = 1:pwEng.numberCashFlowTimes
      pwEng.numberCashFlowsThisIndex[i][j] = 0
      for k = 1:pwEng.numberRates + 1
        pwEng.totalCashFlowsThisIndex[i][j, k] = 0.0
      end
    end
    for l = 1:pwEng.numberRates
      for m = 1:pwEng.numberSteps + 1
        pwEng.V[i][m, l] = 0.0
      end
    end
  end

  weight = start_new_path!(pwEng.evolver)
  reset!(pwEng.product)

  isDone = false
  thisStep = 1
  while ~isDone
    thisStep = pwEng.evolver.currentStep
    storeStep = thisStep + 1
    weight *= advance_step!(pwEng.evolver)

    isDone = next_time_step!(pwEng.product, current_state(pwEng.evolver), pwEng.numberCashFlowsThisStep, pwEng.cashFlowsGenerated)

    pwEng.lastForwards = copy(pwEng.currentForwards)

    pwEng.currentForwards = copy(current_state(pwEng.evolver).forwardRates)

    for i = 1:pwEng.numberRates
      x = discount_ratio(current_state(pwEng.evolver), i+1, i)
      pwEng.stepsDiscounts[i+1] = x
      pwEng.StepsDiscountsSquared[storeStep, i] = x*x

      pwEng.LIBORRatios[storeStep, i] = pwEng.currentForwards[i] / pwEng.lastForwards[i]
      pwEng.LIBORRates[storeStep, i] = pwEng.currentForwards[i]
      pwEng.Discounts[storeStep, i+1] = discount_ratio(current_state(pwEng.evolver), i+1, 1)
    end

    get_bumps!(pwEng.jacobianComputers[thisStep], pwEng.lastForwards, pwEng.stepsDiscounts, pwEng.currentForwards, brownians_this_step(pwEng.evolver),
              pwEng.jacobiansThisPaths[thisStep])

    # for each product and each cash flow
    for i = 1:pwEng.numberProducts, j = 1:pwEng.numberCashFlowsThisStep[i]
      k = pwEng.cashFlowsGenerated[i][j].timeIndex
      pwEng.numberCashFlowsThisIndex[i][k] += 1

      for l = 1:pwEng.numberRates + 1
        pwEng.totalCashFlowsThisIndex[i][k, l] += pwEng.cashFlowsGenerated[i][j].amount[l] * weight
      end
    end

    # if done
    #   break
    # end
  end

  # ok we've gathered cash-flows, still have to do backwards computation
  factors = pwEng.pseudoRootStructure.numberOfFactors
  taus = pwEng.pseudoRootStructure.evolution.rateTaus

  flowsFound = false
  finalStepDone = thisStep

  for currentStep = pwEng.numberSteps:-1:1
    stepToUse = min(currentStep, finalStepDone) + 1
    for k in eachindex(pwEng.cashFlowIndicesThisStep[currentStep])
      cashFlowIndex = pwEng.cashFlowIndicesThisStep[currentStep][k]

      # first check to see if anything actually happened before spending time on computing stuff
      noFlows = true
      for l = 1:pwEng.numberProducts
        noFlows = noFlows && (pwEng.numberCashFlowsThisIndex[l][cashFlowIndex] == 0)
        if ~noFlows
          break
        end
      end

      flowsFound = flowsFound || ~noFlows

      if ~noFlows
        if pwEng.doDeflation
          pwEng.deflatorAndDerivatives = get_factors(pwEng.discounters[cashFlowIndex], pwEng.LIBORRates, pwEng.Discounts, stepToUse)
        end

        for j = 1:pwEng.numberProducts
          if pwEng.numberCashFlowsThisIndex[j][cashFlowIndex] > 0
            deflatedCashFlow = pwEng.totalCashFlowsThisIndex[j][cashFlowIndex, 1]
            if pwEng.doDeflation
              deflatedCashFlow *= pwEng.deflatorAndDerivatives[1]
            end

            pwEng.numerairesHeld[j] += deflatedCashFlow

            for i = 2:pwEng.numberRates + 1
              thisDerivative = pwEng.totalCashFlowsThisIndex[j][cashFlowIndex, i]
              if pwEng.doDeflation
                thisDerivative *= pwEng.deflatorAndDerivatives[1]
                thisDerivative += pwEng.totalCashFlowsThisIndex[j][cashFlowIndex, 1] * pwEng.deflatorAndDerivatives[i]
              else
                pwEng.fullDerivatives[i-1] = thisDerivative
              end

              pwEng.V[j][stepToUse, i-1] += thisDerivative
            end # end of for i = 2:numberRates + 1
          end # end of (numberCashFlowsThisIndex[j][cashFlowIndex] > 0)
        end # end of j = 1:numberProducts
      end # end of ~noFlows
    end # end of k in eachindex(cashFlowIndicesThisStep)

    # need to do backwards updating
    if flowsFound
      nextStepToUse = min(currentStep - 1, finalStepDone)
      nextStepIndex = nextStepToUse + 1

      if nextStepIndex != stepToUse # then we need to update V
        thisPseudoRoot = pwEng.pseudoRootStructure.pseudoRoots[currentStep]

        for i = 1:pwEng.numberProducts
          # compute partials
          for f = 1:factors
            libor = pwEng.LIBORRates[stepToUse, pwEng.numberRates]
            V = pwEng.V[i][stepToUse, pwEng.numberRates]
            pseudo = thisPseudoRoot[pwEng.numberRates, f]
            thisPartialTerm = libor * V * pseudo
            pwEng.partials[f, pwEng.numberRates] = thisPartialTerm

            for r = pwEng.numberRates-1:-1:1
              thisPartialTermr = pwEng.LIBORRates[stepToUse, r] * pwEng.V[i][stepToUse, r] * thisPseudoRoot[r, f]
              pwEng.partials[f, r] = pwEng.partials[f, r+1] + thisPartialTermr
            end
          end # end of f = 1:factors

          for j = 1:pwEng.numberRates
            nextV = pwEng.V[i][stepToUse, j] * pwEng.LIBORRatios[stepToUse, j]
            pwEng.V[i][nextStepIndex, j] = nextV

            summandTerm = 0.0
            for f = 1:factors
              summandTerm += thisPseudoRoot[j, f] * pwEng.partials[f, j]
            end

            summandTerm *= taus[j] * pwEng.StepsDiscountsSquared[stepToUse, j]
            pwEng.V[i][nextStepIndex, j] += summandTerm
          end # end of for j = 1:numberRates
        end # end of i = 1:numberProducts
      end # end of if nextStepIndex != stepToUse
    end # end of if flowsFound
  end # end of for currentStep = numberSteps:-1:1

  # all V matricies computed, we now compute the elementary vegas for this path
  for i = 1:pwEng.numberProducts, j = 1:pwEng.numberSteps
    nextIndex = j+1

    # we know V, we need to pair against the sensitivity of the rate to the elementary vega
    # note the simplification here arising from the fact that the elementary vega affects the evolution on precisely one step
    for k = 1:pwEng.numberRates, f = 1:pwEng.factors
      sensitivity = 0.0
      for r = 1:pwEng.numberRates
        sensitivity += pwEng.V[i][nextIndex, r] * pwEng.jacobiansThisPaths[j][r][k, f]
      end

      pwEng.elementaryVegasThisPath[i][j][k, f]
    end
  end

  # write answer into values
  entriesPerProduct = 1 + pwEng.numberRates + pwEng.numberElementaryVegas

  for i = 1:pwEng.numberProducts
    values[((i - 1) * entriesPerProduct) + 1] = pwEng.numerairesHeld[i] * pwEng.initialNumeraireValue

    for j = 1:pwEng.numberRates
      values[(i - 1) * entriesPerProduct + 1 + j] = pwEng.V[i][1, j] * pwEng.initialNumeraireValue
    end

    for k = 1:pwEng.numberSteps, l = 1:pwEng.numberRates, m = 1:pwEng.factors
      values[(i-1) * entriesPerProduct + pwEng.numberRates + 1 + m + (l - 1) * factors + (k - 1) * pwEng.numberRates * pwEng.factors] =
            pwEng.elementaryVegasThisPath[i][k][l, m] * pwEng.initialNumeraireValue
    end
  end

  return 1.0 # we have put the weight in already, this results in lower variance since weight changes along the path
end

function multiple_path_values_elementary!(pwEng::PathwiseVegasOuterAccountingEngine, means::Vector{Float64}, errors::Vector{Float64}, numberOfPaths::Int)
  numberOfElementaryVegas = pwEng.numberRates * pwEng.numberSteps * pwEng.factors
  values = Vector{Float64}(number_of_products(pwEng.product) * (1 + pwEng.numberRates + numberOfElementaryVegas))
  resize!(means, length(values))
  resize!(errors, length(values))

  sums = zeros(length(values))
  sumsqs = zeros(length(values))

  for i = 1:numberOfPaths
    single_path_values!(pwEng, values)
    for j in eachindex(values)
      sums[j] += values[j]
      sumsqs[j] += values[j] * values[j]
    end
  end

  for j in eachindex(values)
    means[j] = sums[j] / numberOfPaths
    meanSq = sumsqs[j] / numberOfPaths
    variance = meanSq - means[j] * means[j]
    errors[j] = sqrt(variance / numberOfPaths)
  end

  return means, errors
end

function multiple_path_values!(pwEng::PathwiseVegasOuterAccountingEngine, means::Vector{Float64}, errors::Vector{Float64}, numberOfPaths::Int)
  allMeans = Vector{Float64}()
  allErrors = Vector{Float64}()

  multiple_path_values_elementary!(pwEng, allMeans, allErrors, numberOfPaths)

  outDataPerProduct = 1 + pwEng.numberRates + pwEng.numberBumps
  inDataPerProduct = 1 + pwEng.numberRates + pwEng.numberElementaryVegas

  resize!(means, (1 + pwEng.numberRates + pwEng.numberBumps) * pwEng.numberProducts)
  resize!(errors, (1 + pwEng.numberRates + pwEng.numberBumps) * pwEng.numberProducts)

  for p = 1:pwEng.numberProducts, i = 1:pwEng.numberRates + 1
    means[i + (p - 1) * outDataPerProduct] = allMeans[i + (p - 1) * outDataPerProduct]
    errors[i + (p - 1) * outDataPerProduct] = allErrors[i + (p - 1) * outDataPerProduct]

    for bump = 1:pwEng.numberBumps
      thisVega = 0.0
      for t = 1:pwEng.numberSteps, r = 1:pwEng.numberRates, f = 1:pwEng.factors
        thisVega += pwEng.vegaBumps[t][bump][r, f] * allMeans[(p - 1) * inDataPerProduct + 1 + pwEng.numberRates + (t - 1) * pwEng.numberRates * pwEng.factors + (r - 1) * pwEng.factors + f]
      end

      means[(p - 1) * outDataPerProduct + 1 + pwEng.numberRates + bump] = thisVega
    end
  end

  return means, errors
end
