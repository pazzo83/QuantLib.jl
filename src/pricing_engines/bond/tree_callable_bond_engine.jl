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

function _calculate!(pe::TreeCallibleFixedRateEngine, bond::CallableFixedRateBond)
  tsmodel = pe.model

  refDate = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  callableBond = DiscretizedCallableFixedRateBond(bond, refDate, dc)

  if isdefined(pe, :common)
    lattice = pe.common.lattice
  else
    times = mandatory_times(callableBond)
    tg = TimeGrid(times, pe.timeSteps)
    lattice = tree(pe.model, tg)
    pe.common = LatticeShortRateModelEngineCommon(tg, lattice)
  end

  redemptionTime = year_fraction(dc, refDate, callableBond.args.redemptionDate)

  initialize!(callableBond, lattice.treeLattice, redemptionTime)
  rollback!(callableBond, 0.0)
  bond.settlementValue = present_value(callableBond)

  return pe, bond
end
