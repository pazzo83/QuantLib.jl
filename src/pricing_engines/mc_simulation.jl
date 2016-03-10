type MCSimulation{RSG <: AbstractRandomSequenceGenerator, T <: MCTrait}
  antitheticVariate::Bool
  controlVariate::Bool
  rsg::RSG
  mcTrait::T
  mcModel::MonteCarloModel

  MCSimulation{RSG, T}(antitheticVariate::Bool, controlVariate::Bool, rsg::RSG, mcTrait::T) =
              new{RSG, T}(antitheticVariate, controlVariate, rsg, mcTrait)
end

max_error(err::Float64) = err

function get_value!(mcsim::MCSimulation, tolerance::Float64, maxSamples::Int = typemax(Int), minSamples::Int = 1023)
  sampleNumber = size(mcsim.mcModel.sampleAccumulator.samplesMatrix)[1]
  if sampleNumber < minSamples
    idxNum = sampleNumber + 1
    add_samples!(mcsim.mcModel, minSamples - sampleNumber, idxNum) # add an index so we can add to this
    sampleNumber = size(mcsim.mcModel.sampleAccumulator.samplesMatrix)[1]
  end

  err = error_estimate(mcsim.mcModel.sampleAccumulator)
  while max_error(err) > tolerance
    sampleNumber < maxSamples || error("max number of samples reached while error still above tolerance")
    # conservative estimate of how many samples are needed
    order = max_error(err * err) / tolerance / tolerance
    nextBatch = round(Int, floor(max(sampleNumber * order * 0.8 - sampleNumber, minSamples)))
    nextBatch = min(nextBatch, maxSamples - sampleNumber)
    idxNum = sampleNumber + 1
    sampleNumber += nextBatch
    add_samples!(mcsim.mcModel, nextBatch, idxNum)
    err = error_estimate(mcsim.mcModel.sampleAccumulator)
  end
  return stats_mean(mcsim.mcModel.sampleAccumulator)
end

function get_value_with_samples!(mcsim::MCSimulation, samples::Int)
  sampleNumber = size(mcsim.mcModel.sampleAccumulator.samplesMatrix)[1]

  samples >= sampleNumber || error("number of already simulated samples greater than requested samples")

  idxNum = sampleNumber + 1
  add_samples!(mcsim.mcModel, samples - sampleNumber, idxNum)

  return stats_mean(mcsim.mcModel.sampleAccumulator)
end

function _calculate!(mcsim::MCSimulation, pe::PricingEngine, inst::Instrument, requiredTolerance::Float64, requiredSamples::Int, maxSamples::Int)
  # TODO check if control variate
  mcsim.mcModel = MonteCarloModel(path_generator(pe, inst), path_pricer(pe, inst), gen_RiskStatistics(), mcsim.antitheticVariate)

  if requiredTolerance != -1.0
    if maxSamples != -1
      get_value!(mcsim, requiredTolerance, maxSamples)
    else
      get_value!(mcsim, requiredTolerance)
    end
  else
    get_value_with_samples!(mcsim, requiredSamples)
  end

  return mcsim
end
