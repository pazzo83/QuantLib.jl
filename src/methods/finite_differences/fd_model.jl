## Main Finite Difference Model ##
type FiniteDifferenceModel{T}
  evolver::T
  stoppingTimes::Vector{Float64}

  FiniteDifferenceModel{T}(evolver::T, stoppingTimes::Vector{Float64}) = new(evolver, sort(unique(stoppingTimes)))
end

function rollback_impl!(model::FiniteDifferenceModel, a::Vector{Float64}, from::Float64, to::Float64, steps::Int, condition::StepCondition)
  dt = (from - to) / steps
  t = from
  set_step!(model.evolver, dt)

  if length(model.stoppingTimes) > 0 && model.stoppingTimes[end] == from
    apply_to!(condition, a, from)
  end

  for i = 1:steps
    _now = t
    _next = t - dt
    hit = false
    for j = length(model.stoppingTimes):-1:1
      if _next <= model.stoppingTimes[j] && model.stoppingTimes[j] < _now
        # a stopping time was hit
        hit = true

        # perform a small step to stoppingTimes[j]
        set_step!(model.evolver, _now - model.stoppingTimes[j])
        step!(model.evolver, a, _now)

        apply_to!(condition, a, model.stoppingTimes[j])

        # and continue the cycle
        _now = model.stoppingTimes[j]
      end
    end

    # if we did hit
    if hit
      # we might have to make a small step to complete the big one
      if _now > _next
        set_step!(model.evolver, _now - _next)
        step!(model.evolver, a, _now)
        apply_to!(condition, a, _next)
      end

      # and in any case, we have to reset the evolver to the default step
      set_step!(model.evolver, dt)
    else
      # if we didn't, the evolver is already set to the default step, which is ok
      step!(model.evolver, a, _now)
      apply_to!(condition, a, _next)
    end
    t -= dt
  end

  return model, a
end

rollback!(model::FiniteDifferenceModel, a::Vector{Float64}, from::Float64, to::Float64, steps::Int, condition::StepCondition) = rollback_impl!(model, a, from, to, steps, condition)
