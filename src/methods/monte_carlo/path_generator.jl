type PathGenerator{RSG <: AbstractRandomSequenceGenerator, I <: Integer, S <: StochasticProcess1D}
  brownianBridge::Bool
  generator::RSG
  dimension::I
  timeGrid::TimeGrid
  process::S
  nextSample::Sample
  temp::Float64
  bb::BrownianBridge
end

function PathGenerator(process::StochasticProcess, tg::TimeGrid, generator::AbstractRandomSequenceGenerator, brownianBridge::Bool)
  generator.dimension == length(tg.times) - 1 || error("wrong dimensions")

  return PathGenerator(brownianBridge, generator, generator.dimension, tg, process, Sample(Path(tg), 1.0), zeros(dimension), BrownianBridge(tg))
end

get_next!(pg::PathGenerator) = get_next!(pg, false)

function get_next!(pg::PathGenerator, antithetic::Bool)
  sequenceVals, sequenceWeight = antithetic ? last_sequence(generator) : next_sequence!(generator)

  # TODO Brownian Bridge
  temp = copy(sequenceVals)

  nextSample.weight = sequenceWeight
  path = next.value # TODO check to see if this is just a pass by reference

  path[1] = get_x0(pg.process)

  for i = 2:length(path)
    t = pg.timeGrid[i-1]
    dt = pg.timeGrid.dt[i - 1]
    path[i] = evolve(process, path[i-1], dt, antithetic ? -temp[i-1] : temp[i-1])
  end

  return next.value
end
