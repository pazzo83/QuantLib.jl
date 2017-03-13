type SubProduct{M <: MarketModelMultiProduct}
  product::M
  multiplier::Float64
  numberOfCashflows::Vector{Int}
  cashFlows::Vector{Vector{MarketModelCashFlow}}
  timeIndices::Vector{Int}
  isDone::Bool
end

# product is already cloned
SubProduct(product::MarketModelMultiProduct, multiplier::Float64, isDone::Bool) =
          SubProduct(product, multiplier, Vector{Int}(), Vector{Vector{MarketModelCashFlow}}(), Vector{Int}(), isDone)


type MarketModelComposite <: MarketModelMultiProduct
  components::Vector{SubProduct}
  rateTimes::Vector{Float64}
  evolutionTimes::Vector{Float64}
  finalized::Bool
  currentIndex::Int
  cashFlowTimes::Vector{Float64}
  allEvolutionTimes::Vector{Vector{Float64}}
  isInSubset::Vector{BitArray{1}}
  evolution::EvolutionDescription
end

MarketModelComposite() = MarketModelComposite(Vector{SubProduct}(), Vector{Float64}(), Vector{Float64}(), false, 1, Vector{Float64}(), Vector{Vector{Float64}}(),
                            Vector{BitArray{1}}(), EvolutionDescription())

function add_product!(mm::MarketModelComposite, product::MarketModelMultiProduct, multiplier::Float64 = 1.0)
  product = clone(product)
  d = get_evolution(product)
  if ~isempty(mm.components)
    # enforce pre-conditions
    d1 = get_evolution(mm.components[1].product)
    rateTimes1 = d1.rateTimes
    rateTimes2 = d.rateTimes
    (length(rateTimes1) == length(rateTimes2) && rateTimes1 == rateTimes2) || error("incompatible rate times")
  end

  push!(mm.components, SubProduct(product, multiplier, false))
  push!(mm.allEvolutionTimes, d.evolutionTimes)

  return mm
end

function finalize!(mm::MarketModelComposite)
  ~mm.finalized || error("already finalized")
  ~isempty(mm.components) || error("no subproducts provided")

  # fetch the rate times from the first subproduct (we already checked that they are the same)
  description = get_evolution(mm.components[1].product)
  mm.rateTimes = description.rateTimes

  mm.evolutionTimes, mm.isInSubset = merge_times(mm.allEvolutionTimes)

  allCashFlowTimes = Vector{Float64}()

  # now for each subproduct
  @inbounds @simd for i in eachindex(mm.components)
    d = get_evolution(mm.components[i].product)

    # collect all possible cash flow times
    cashFlowTimes = possible_cash_flow_times(mm.components[i].product)
    append!(allCashFlowTimes, cashFlowTimes)

    # allocate working vectors
    n = number_of_products(mm.components[i].product)
    mm.components[i].numberOfCashflows = Vector{Int}(n)
    mm.components[i].cashFlows = Vector{MarketModelCashFlow}[
                                [MarketModelCashFlow() for i=1:max_number_of_cashflows_per_product_per_step(mm.components[i].product)] for j = 1:n]
  end

  # all information collected, we can sort and compact the vector of all cashflow times
  sort!(allCashFlowTimes)
  allCashFlowTimes = unique(allCashFlowTimes)

  append!(mm.cashFlowTimes, allCashFlowTimes)

  # and map each product's cash flow time into the total vector
  @inbounds @simd for i in eachindex(mm.components)
    productTimes = possible_cash_flow_times(mm.components[i].product)
    mm.components[i].timeIndices = Int[findfirst(mm.cashFlowTimes, productTimes[j]) for j in eachindex(productTimes)]
  end

  mm.evolution = EvolutionDescription(mm.rateTimes, mm.evolutionTimes)

  mm.finalized = true

  return mm
end

function next_time_step!(mm::MarketModelComposite, currentState::CurveState, numberCashFlowsThisStep::Vector{Int}, cashFlowsGenerated::Vector{Vector{MarketModelCashFlow}})
  mm.finalized || error("composite not finalized")

  isDone = true
  n = 1
  offset = 0
  @inbounds @simd for i in eachindex(mm.components)
    if mm.isInSubset[n][mm.currentIndex] && ~mm.components[i].isDone
      # make it evolve
      thisDone = next_time_step!(mm.components[i].product, currentState, mm.components[i].numberOfCashflows, mm.components[i].cashFlows)

      # and copy the results. Time indices need to be remapped so that they point into all cash-flow times
      # Amounts need to be adjusted by the corresponding multiplier
      for j = 1:number_of_products(mm.components[i].product)
        numberCashFlowsThisStep[j+offset] = mm.components[i].numberOfCashflows[j]
        for k = 1:mm.components[i].numberOfCashflows[j]
          cfFrom = mm.components[i].cashFlows[j][k]
          cfTo = cashFlowsGenerated[j+offset][k]
          cfTo.timeIndex = mm.components[i].timeIndices[cfFrom.timeIndex]
          cfTo.amount = cfFrom.amount * mm.components[i].multiplier
        end
      end

      # finally set isDone to false if this product isn't done
      isDone = isDone && thisDone
    else
      for j = 1:number_of_products(mm.components[i].product)
        numberCashFlowsThisStep[j+offset] = 0
      end
    end

    # the offset is updated whether or not the product was evolved
    offset += number_of_products(mm.components[i].product)
    n += 1
  end
  mm.currentIndex += 1
  return isDone
end


function number_of_products(mm::MarketModelComposite)
  res = 0
  @simd for i in eachindex(mm.components)
    @inbounds res += number_of_products(mm.components[i].product)
  end

  return res
end

function max_number_of_cashflows_per_product_per_step(mm::MarketModelComposite)
  res = 0
  @simd for i in eachindex(mm.components)
    @inbounds res  = max(res, max_number_of_cashflows_per_product_per_step(mm.components[i].product))
  end

  return res
end

function possible_cash_flow_times(mm::MarketModelComposite)
  mm.finalized || error("composite not finalized")
  return mm.cashFlowTimes
end

function get_evolution(mm::MarketModelComposite)
  mm.finalized || error("composite not finalized")
  return mm.evolution
end

function reset!(mm::MarketModelComposite)
  @inbounds @simd for i in eachindex(mm.components)
    reset!(mm.components[i].product)
    mm.components[i].isDone = false
  end
  mm.currentIndex = 1
  return mm
end

clone(mm::MarketModelComposite) = MarketModelComposite(deepcopy(mm.components), copy(mm.rateTimes), copy(mm.evolutionTimes), mm.finalized, mm.currentIndex,
                                    copy(mm.cashFlowTimes), deepcopy(mm.allEvolutionTimes), deepcopy(mm.isInSubset), clone(mm.evolution))
