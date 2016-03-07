type MCSimulation{RNG <: AbstractRandomNumberGenerator, T <: MCTrait}
  antitheticVariate::Bool
  controlVariate::Bool
  rng::RNG
  stats::Statistics
  mcTrait::T
  mcModel::MonteCarloModel

  MCSimulation(antitheticVariate::Bool, controlVariate::Bool, rng::RNG, stats::Statistics, mcTrait::T) =
              new(antitheticVariate, controlVariate, rng, stats, mcTrait)
end
