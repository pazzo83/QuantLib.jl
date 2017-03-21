type TreeCallableFixedRateEngine{S <: ShortRateModel, T <: ShortRateTree} <: LatticeShortRateModelEngine{S}
  model::S
  timeSteps::Int
  common::LatticeShortRateModelEngineCommon{T}
  latticeBuilt::Bool

  function TreeCallableFixedRateEngine{S, T}(model::S, timeSteps::Int, common::LatticeShortRateModelEngineCommon{T}, latticeGen::Bool = true)
    te = new{S, T}(model, timeSteps, common, latticeGen)

    add_observer!(model, te)

    return te
  end
end

function TreeCallableFixedRateEngine{S <: ShortRateModel}(model::S, timeSteps::Int)
  # create empty tree
  tg = TimeGrid([1.0], 1)
  lattice = tree(model, tg)
  common = LatticeShortRateModelEngineCommon{typeof(lattice)}(tg, lattice)
  TreeCallableFixedRateEngine{S, typeof(lattice)}(model, timeSteps, common, false)
end

function _calculate!(pe::TreeCallableFixedRateEngine, bond::CallableFixedRateBond)
  tsmodel = pe.model

  refDate = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  callableBond = DiscretizedCallableFixedRateBond(bond, refDate, dc, pe.common.lattice.treeLattice)

  if ~pe.latticeBuilt
    times = mandatory_times(callableBond)
    tg = TimeGrid(times, pe.timeSteps)
    pe.common.tg = tg
    update!(pe)
    pe.latticeBuilt = true
  end

  # if isdefined(pe, :common)
  #   lattice = pe.common.lattice
  # else
  #   times = mandatory_times(callableBond)
  #   tg = TimeGrid(times, pe.timeSteps)
  #   lattice = tree(pe.model, tg)
  #   pe.common = LatticeShortRateModelEngineCommon(tg, lattice)
  # end

  lattice = pe.common.lattice

  redemptionTime = year_fraction(dc, refDate, callableBond.args.redemptionDate)

  initialize!(callableBond, lattice.treeLattice, redemptionTime)
  rollback!(callableBond, 0.0)
  bond.settlementValue = present_value(callableBond)

  return pe, bond
end
