type TreeCallibleFixedRateEngine{S <: ShortRateModel} <: LatticeShortRateModelEngine{S}
  model::S
  timeSteps::Int
  common::LatticeShortRateModelEngineCommon

  function TreeCallibleFixedRateEngine{S}(model::S, timeSteps::Int)
    te = new{S}(model, timeSteps)

    add_observer!(model, te)

    return te
  end
end

TreeCallibleFixedRateEngine{S <: ShortRateModel}(model::S, timeSteps::Int) = TreeCallibleFixedRateEngine{S}(model, timeSteps)
