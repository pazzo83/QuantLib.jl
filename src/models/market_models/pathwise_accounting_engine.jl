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
    idx = upper_bound(evolutionTimes, cashFlowTimes[i])
    if idx != 1
      idx -= 1
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
