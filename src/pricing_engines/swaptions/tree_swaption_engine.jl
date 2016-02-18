type LatticeShortRateModelEngineCommon{T <: ShortRateTree}
  tg::TimeGrid
  lattice::T
end


type TreeSwaptionEngine{S <: ShortRateModel, I <: Integer, Y <: YieldTermStructure} <: LatticeShortRateModelEngine{S, Y}
  model::S
  timeSteps::I
  common::LatticeShortRateModelEngineCommon
  ts::Y

  function call{S, I}(::Type{TreeSwaptionEngine}, m::S, tsteps::I)
    t = new{S, I, YieldTermStructure}(m, tsteps)
    add_observer!(m, t)

    return t
  end

  call{S, I}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, l::LatticeShortRateModelEngineCommon) = new{S, I, YieldTermStructure}(m, tsteps, l)

  call{S, I, Y}(::Type{TreeSwaptionEngine}, m::S, tsteps::I, l::LatticeShortRateModelEngineCommon, ts::Y) = new{S, I, T, Y}(m, tsteps, l, ts)
end

function TreeSwaptionEngine{S <: ShortRateModel}(model::S, tg::TimeGrid)
  lattice = tree(model, tg)
  ts = TreeSwaptionEngine(model, 0, LatticeShortRateModelEngineCommon(tg, lattice))

  add_observer!(model, ts)

  return ts
end

function update!(eng::LatticeShortRateModelEngine)
  if length(eng.common.tg.times) > 0
    eng.common.lattice = tree(eng.model, eng.common.tg)
  end

  return eng
end

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
  @simd for i = 1:length(stoppingTimes)
    @inbounds stoppingTimes[i] = year_fraction(dc, refDate, swaption.exercise.dates[i])
  end

  initialize!(dSwaption, lattice.treeLattice, stoppingTimes[end])

  nextExerciseIdx = findnext(greater_than_or_equal_to, stoppingTimes, 1, 0.0)

  nextExercise = nextExerciseIdx != 0 ? stoppingTimes[nextExerciseIdx] : stoppingTimes[end]

  rollback!(dSwaption, nextExercise)
  swaption.results.value = present_value(dSwaption)
end
