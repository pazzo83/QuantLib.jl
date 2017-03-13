type AccountingEngine{E <: AbstractMarketModelEvolver, P <: MarketModelMultiProduct}
  evolver::E
  product::P
  initialNumeraireValue::Float64
  numberProducts::Int
  numerairesHeld::Vector{Float64}
  numberCashFlowsThisStep::Vector{Int}
  cashFlowsGenerated::Vector{Vector{MarketModelCashFlow}}
  discounters::Vector{MarketModelDiscounter}
end

function AccountingEngine(evolver::AbstractMarketModelEvolver, product::MarketModelMultiProduct, initialNumeraireValue::Float64)
  product = clone(product)
  numProd = number_of_products(product)

  cashFlowsGenerated = Vector{MarketModelCashFlow}[[MarketModelCashFlow() for i = 1:max_number_of_cashflows_per_product_per_step(product)] for j = 1:numProd]

  cashFlowTimes = possible_cash_flow_times(product)
  rateTimes = get_evolution(product).rateTimes

  discounters = MarketModelDiscounter[MarketModelDiscounter(cashFlowTimes[i], rateTimes) for i in eachindex(cashFlowTimes)]

  return AccountingEngine(evolver, product, initialNumeraireValue, numProd, Vector{Float64}(numProd), Vector{Int}(numProd), cashFlowsGenerated, discounters)
end

function single_path_values!(ae::AccountingEngine, vals::Vector{Float64})
  fill!(ae.numerairesHeld, 0.0)
  weight = start_new_path!(ae.evolver)
  reset!(ae.product)

  principalInNumerairePortfolio = 1.0

  isDone = false

  while true
    thisStep = ae.evolver.currentStep
    weight *= advance_step!(ae.evolver)
    # println("forwards: ", ae.evolver.forwards)
    # println("logs: ", ae.evolver.logForwards)
    isDone = next_time_step!(ae.product, current_state(ae.evolver), ae.numberCashFlowsThisStep, ae.cashFlowsGenerated)

    numeraire = ae.evolver.numeraires[thisStep]

    # println("cfs: ", ae.cashFlowsGenerated)

    # for each product
    @inbounds @simd for i = 1:ae.numberProducts
      # and each cash flow
      cfs = ae.cashFlowsGenerated[i]

      for j = 1:ae.numberCashFlowsThisStep[i]
        # convert the cash flows to numeraires
        # this is done by calculating the number of numeraire bonds corresponding to such cash flow
        discounter = ae.discounters[cfs[j].timeIndex]
        bonds = cfs[j].amount * numeraire_bonds(discounter, current_state(ae.evolver), numeraire)

        # and adding the newly bought bonds to the number of numeraires held
        ae.numerairesHeld[i] += bonds / principalInNumerairePortfolio
      end
    end

    # error("die")

    if ~isDone
      # The numeraire might change between steps. This implies
      # that we might have to convert the numeraire bonds for
      # this step into a corresponding amount of numeraire
      # bonds for the next step. This can be done by changing
      # the principal of the numeraire and updating the number
      # of bonds in the numeraire portfolio accordingly.
      nextNumeraire = ae.evolver.numeraires[thisStep + 1]
      principalInNumerairePortfolio *= discount_ratio(current_state(ae.evolver), numeraire, nextNumeraire)
    end

    if isDone
      break
    end
  end

  @simd for i in eachindex(ae.numerairesHeld)
    @inbounds vals[i] = ae.numerairesHeld[i] * ae.initialNumeraireValue
  end

  return weight
end

function multiple_path_values!(ae::AccountingEngine, stats::GenericSequenceStats, numberOfPaths::Int)
  vals = Vector{Float64}(number_of_products(ae.product))
  QuantLib.Math.reset!(stats, length(vals), numberOfPaths)

  @inbounds @simd for i = 1:numberOfPaths
    weight = single_path_values!(ae, vals)
    # println("vals ", vals)
    # if i > 2
    #   error("DIE")
    # end
    add_sample!(stats, vals, i, weight)
  end

  return stats
end
