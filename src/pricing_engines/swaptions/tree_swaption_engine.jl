type TreeSwaptionEngine{S <: ShortRateModel, T <: ShortRateTree} <: LatticeShortRateModelEngine{S}
  model::S
  timeSteps::Int
  common::LatticeShortRateModelEngineCommon{T}
  latticeBuilt::Bool
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
  function TreeSwaptionEngine{S, T}(model::S, timeSteps::Int, common::LatticeShortRateModelEngineCommon{T}, latticeGen::Bool = true)
    ts = new{S, T}(model, timeSteps, common, latticeGen)
    add_observer!(model, ts)

    return ts
  end

  # function TreeSwaptionEngine{S}(model::S, timeSteps::Int)
  #   # create empty tree
  #   tg = TimeGrid([1.0], 1)
  #   lattice = tree(model, tg)
  #   common = LatticeShortRateModelEngineCommon{typeof(lattice)}(tg, lattice)
  #   ts = new{S, typeof(lattice)}(model, timeSteps, common, false)
  #   add_observer!(model, ts)
  #
  #   return ts
  # end
end

function TreeSwaptionEngine{S <: ShortRateModel}(model::S, tg::TimeGrid)
  lattice = tree(model, tg)
  return TreeSwaptionEngine{S, typeof(lattice)}(model, 0, LatticeShortRateModelEngineCommon(tg, lattice))
end

function TreeSwaptionEngine{S <: ShortRateModel}(model::S, timeSteps::Int)
  # create empty tree
  tg = TimeGrid([1.0], 1)
  lattice = tree(model, tg)
  common = LatticeShortRateModelEngineCommon{typeof(lattice)}(tg, lattice)
  return TreeSwaptionEngine{S, typeof(lattice)}(model, timeSteps, common, false)
end

# methods
function _calculate!(pe::TreeSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model

  refDate = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  dSwaption = DiscretizedSwaption(swaption, refDate, dc, pe.common.lattice.treeLattice)

  if ~pe.latticeBuilt
    times = mandatory_times(dSwaption)
    tg = TimeGrid(times, pe.timeSteps)
    pe.common.tg = tg
    update!(pe)
    pe.latticeBuilt = true
  end
  lattice = pe.common.lattice

  stoppingTimes = zeros(length(swaption.exercise.dates))
  # @simd for i in eachindex(stoppingTimes)
  #   @inbounds stoppingTimes[i] = year_fraction(dc, refDate, swaption.exercise.dates[i])
  # end
  map!(x -> year_fraction(dc, refDate, x), stoppingTimes, swaption.exercise.dates)

  initialize!(dSwaption, lattice.treeLattice, stoppingTimes[end])

  nextExerciseIdx = findnext(greater_than_or_equal_to, stoppingTimes, 1, 0.0)

  nextExercise = nextExerciseIdx != 0 ? stoppingTimes[nextExerciseIdx] : stoppingTimes[end]

  rollback!(dSwaption, nextExercise)
  swaption.results.value = present_value(dSwaption)
end
