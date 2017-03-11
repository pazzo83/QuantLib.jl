type TreeSwaptionEngine{S <: ShortRateModel} <: LatticeShortRateModelEngine{S}
  model::S
  timeSteps::Int
  common::LatticeShortRateModelEngineCommon
  # ts::Y

  # function call{S, I}(::Type{TreeSwaptionEngine}, m::S, tsteps::I)
  #   t = new{S, I, YieldTermStructure}(m, tsteps)
  #   add_observer!(m, t)
  #
  #   return t
  # end
  #
  # call{S, I}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, l::LatticeShortRateModelEngineCommon) = new{S, I, YieldTermStructure}(m, tsteps, l)
  #
  # call{S, I, Y}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, l::LatticeShortRateModelEngineCommon, ts::Y) = new{S, I, T, Y}(m, tsteps, l, ts)
  function TreeSwaptionEngine{S}(model::S, timeSteps::Int, common::LatticeShortRateModelEngineCommon)
    ts = new{S}(model, timeSteps, common)
    add_observer!(model, ts)

    return ts
  end

  function TreeSwaptionEngine{S}(model::S, timeSteps::Int)
    ts = new{S}(model, timeSteps)
    add_observer!(model, ts)

    return ts
  end
end

function TreeSwaptionEngine{S <: ShortRateModel}(model::S, tg::TimeGrid)
  lattice = tree(model, tg)
  return TreeSwaptionEngine{S}(model, 0, LatticeShortRateModelEngineCommon(tg, lattice))
end

TreeSwaptionEngine{S <: ShortRateModel}(model::S, timeSteps::Int) = TreeSwaptionEngine{S}(model, timeSteps)

# methods
function _calculate!(pe::TreeSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model

  refDate = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  dSwaption = DiscretizedSwaption(swaption, refDate, dc)

  if isdefined(pe, :common)
    lattice = pe.common.lattice
  else
    times = mandatory_times(dSwaption)
    tg = TimeGrid(times, pe.timeSteps)
    lattice = tree(pe.model, tg)
    pe.common = LatticeShortRateModelEngineCommon(tg, lattice)
  end

  stoppingTimes = zeros(length(swaption.exercise.dates))
  @simd for i in eachindex(stoppingTimes)
    @inbounds stoppingTimes[i] = year_fraction(dc, refDate, swaption.exercise.dates[i])
  end

  initialize!(dSwaption, lattice.treeLattice, stoppingTimes[end])

  nextExerciseIdx = findnext(greater_than_or_equal_to, stoppingTimes, 1, 0.0)

  nextExercise = nextExerciseIdx != 0 ? stoppingTimes[nextExerciseIdx] : stoppingTimes[end]

  rollback!(dSwaption, nextExercise)
  swaption.results.value = present_value(dSwaption)
end
