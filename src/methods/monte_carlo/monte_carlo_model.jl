type SingleVariate <: MCTrait end

type MonteCarloModel{P <: AbstractPathPricer} <: AbstractMonteCarloModel
  pathGenerator::PathGenerator
  pathPricer::P
  sampleAccumulator::RiskStatistics
  isAntitheticVariate::Bool
end

function add_samples!(mcmodel::MonteCarloModel, samples::Int, idx::Int)
  # re-init the risk data
  adding_data!(mcmodel.sampleAccumulator, samples)
  for j = 1:samples
    path = get_next!(mcmodel.pathGenerator)
    price = mcmodel.pathPricer(path.value)

    # TODO Control Variate
    # TODO antithetic variate

    add_sample!(mcmodel.sampleAccumulator, price, path.weight, idx)
    idx += 1
  end

  return mcmodel
end
