struct SingleVariate <: MCTrait end

mutable struct MonteCarloModel{P <: AbstractPathPricer, RSG <: AbstractRandomSequenceGenerator, S <: StochasticProcess1D} <: AbstractMonteCarloModel
  pathGenerator::PathGenerator{RSG, S}
  pathPricer::P
  sampleAccumulator::RiskStatistics
  isAntitheticVariate::Bool
end

function add_samples!(mcmodel::MonteCarloModel, samples::Int, idx::Int=1)
  # re-init the risk data
  adding_data!(mcmodel.sampleAccumulator, samples)
  @inbounds @simd for j = 1:samples
    path = get_next!(mcmodel.pathGenerator)
    price = mcmodel.pathPricer(path.value)

    # TODO Control Variate
    if mcmodel.isAntitheticVariate
      path2 = get_antithetic!(mcmodel.pathGenerator)
      price2 = mcmodel.pathPricer(path2.value)
      add_sample!(mcmodel.sampleAccumulator, (price+price2)/2.0, path.weight, idx)
    else
      add_sample!(mcmodel.sampleAccumulator, price, path.weight, idx)
    end
    idx += 1
  end

  return mcmodel
end
