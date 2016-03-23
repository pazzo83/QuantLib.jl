type MarketModelPathwiseCashRebate <: MarketModelPathwiseMultiProduct
  evolution::EvolutionDescription
  paymentTimes::Vector{Float64}
  amounts::Matrix{Float64}
  numberOfProducts::Int
  currentIndex::Int
end

function MarketModelPathwiseCashRebate(evolution::EvolutionDescription,
                                      paymentTimes::Vector{Float64},
                                      amounts::Matrix{Float64},
                                      numberOfProducts::Int)

  check_increasing_times(paymentTimes)
  amounts_rows, amounts_cols = size(amounts)

  amounts_rows == numberOfProducts || error("the number of rows in the matrix must equal the number of products")
  amounts_cols == length(paymentTimes) || error("the number of cols in the matrix must equal the number of payment times")
  length(evolution.evolutionTimes) == length(paymentTimes) || error("the number of evolution times must equal the number of payment times")

  return MarketModelPathwiseCashRebate(evolution, paymentTimes, amounts, numberOfProducts, 1)
end

possible_cash_flow_times(mm::MarketModelPathwiseCashRebate) = mm.paymentTimes

max_number_of_cashflows_per_product_per_step(::MarketModelPathwiseCashRebate) = 1
