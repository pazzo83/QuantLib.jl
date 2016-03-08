type MCSimulation{RSG <: AbstractRandomSequenceGenerator, T <: MCTrait}
  antitheticVariate::Bool
  controlVariate::Bool
  rsg::RSG
  stats::RiskStatistics
  mcTrait::T
  mcModel::MonteCarloModel

  MCSimulation{RSG, T}(antitheticVariate::Bool, controlVariate::Bool, rsg::RSG, stats::RiskStatistics, mcTrait::T) =
              new{RSG, T}(antitheticVariate, controlVariate, rsg, stats, mcTrait)
end

function get_value!(mcsim::MCSimulation, tolerance::Float64, maxSamples::Int = typemax(Int), minSamples::Int = 1023)
  sampleNumber = size(mcsim.mcModel.sampleAccumulator.samplesMatrix)[1]
  if sampleNumber < minSamples
    idxNum = sampleNumber + 1
    add_samples!(mcsim.mcModel, minSamples - sampleNumber, idxNum) # add an index so we can add to this
    sampleNumber = size(mcsim.mcModel.sampleAccumulator.samplesMatrix)[1]
  end
end

function _calculate!(mcsim::MCSimulation, pe::PricingEngine, inst::Instrument, requiredTolerance::Float64, requiredSamples::Int, maxSamples::Int)
  # TODO check if control variate
  mcsim.mcModel = MonteCarloModel(path_generator(pe, inst), path_pricer(pe, inst), mcsim.stats, mcsim.antitheticVariate)

  if requiredTolerance != -1.0
    if maxSamples != -1
      get_value!(mcsim, requiredTolerance, maxSamples)
    else
      get_value!(mcsim, requiredTolerance)
    end
  end

  return mcsim
end
