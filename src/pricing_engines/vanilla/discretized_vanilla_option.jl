type DiscretizedVanillaOption <: DiscretizedAsset
  stoppingTimes::Vector{Float64}
  args::VanillaOptionArgs
  common::DiscretizedAssetCommon
end

function DiscretizedVanillaOption(vanillaOption::VanillaOption, process::StochasticProcess, grid::TimeGrid)
  stoppingTimes = zeros(length(vanillaOption.exercise.dates))
  for i in eachindex(stoppingTimes)
    stoppingTimes[i] = get_time(process, vanillaOption.exercise.dates[i])

    if ~is_empty(grid)
      # stuff
      stoppingTimes[i] = closest_time(grid, stoppingTimes[i])
    end
  end

  return DiscretizedVanillaOption(stoppingTimes, VanillaOptionArgs(vanillaOption.payoff, vanillaOption.exercise), DiscretizedAssetCommon())
end

function reset!(dvo::DiscretizedVanillaOption, sz::Int)
  dvo.common.values = zeros(sz)
  adjust_values!(dvo)

  return dvo
end

post_adjust_values_impl!(dvo::DiscretizedVanillaOption) = post_adjust_values_impl!(dvo, dvo.args.exercise, dvo.common.time)

function post_adjust_values_impl!(dvo::DiscretizedVanillaOption, ::AmericanExercise, time_now::Float64)
  if time_now <= dvo.stoppingTimes[2] && time_now >= dvo.stoppingTimes[1]
    apply_specific_condition!(dvo)
  end
  return dvo
end

function post_adjust_values_impl!(dvo::DiscretizedVanillaOption, ::EuropeanExercise, time_now::Float64)
  if is_on_time(dvo, dvo.stoppingTimes[1])
    apply_specific_condition!(dvo)
  end
  return dvo
end

function post_adjust_values_impl!(dvo::DiscretizedVanillaOption, ::BermudanExercise, time_now::Float64)
  for i in eachindex(dvo.stoppingTimes)
    if is_on_time(dvo, dvo.stoppingTimes[i])
      apply_specific_condition!(dvo)
    end
  end
  return dvo
end

function apply_specific_condition!(dvo::DiscretizedVanillaOption)
  grid = get_grid(dvo.common.method, dvo.common.time)
  for j in eachindex(dvo.common.values)
    dvo.common.values[j] = max(dvo.common.values[j], dvo.args.payoff(grid[j]))
  end

  return dvo
end
