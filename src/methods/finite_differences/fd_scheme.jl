type Hundsdorfer <: FdmSchemeDescType end
type Douglas <: FdmSchemeDescType end

type FdmSchemeDesc{F <: FdmSchemeDescType}
  schemeType::F
  theta::Float64
  mu::Float64
end

FdmSchemeDesc(t::Hundsdorfer) = FdmSchemeDesc(t, 0.5 + sqrt(3.0) / 6.0, 0.5)
FdmSchemeDesc(t::Douglas) = FdmSchemeDesc(t, 0.5, 0.0)

## Schemes ##
type HundsdorferScheme <: FdScheme
  theta::Float64
  mu::Float64
  map::FdmLinearOpComposite
  bcSet::FdmBoundaryConditionSet
  dt::Float64
end

HundsdorferScheme(theta::Float64, mu::Float64, map::FdmLinearOpComposite, bcSet::FdmBoundaryConditionSet) = HundsdorferScheme(theta, mu, map, bcSet, 0.0)

type DouglasScheme <: FdScheme
  theta::Float64
  map::FdmLinearOpComposite
  bcSet::FdmBoundaryConditionSet
  dt::Float64
end

DouglasScheme(theta::Float64, map::FdmLinearOpComposite, bcSet::FdmBoundaryConditionSet) = DouglasScheme(theta, map, bcSet, 0.0)

set_step!(evolver::FdScheme, dt::Float64) = evolver.dt = dt

function step!(evolver::HundsdorferScheme, a::Vector{Float64}, t::Float64)
  t - evolver.dt > -1e-8 || error("a step towards negative time given")

  set_time!(evolver.map, max(0.0, t - evolver.dt), t)
  set_time!(evolver.bcSet, max(0.0, t - evolver.dt))

  apply_before_applying!(evolver.bcSet, evolver.map)
  y = a + evolver.dt * apply(evolver.map, a)
  apply_after_applying!(evolver.bcSet, y)

  y0 = copy(y)

  for i = 1:get_size(evolver.map)
    rhs = y - evolver.theta * evolver.dt * apply_direction(evolver.map, i, a)
    y = solve_splitting(evolver.map, i, rhs, -evolver.theta * evolver.dt)
  end

  apply_before_applying!(evolver.bcSet, evolver.map)
  yt = y0 + evolver.mu * evolver.dt * apply(evolver.map, y - a)
  apply_after_applying!(evolver.bcSet, yt)

  for i = 1:get_size(evolver.map)
    rhs = yt - evolver.theta * evolver.dt * apply_direction(evolver.map, i, y)
    yt = solve_splitting(evolver.map, i, rhs, -evolver.theta * evolver.dt)
  end

  apply_after_solving!(evolver.bcSet, yt)

  a[:] = yt

  return a, evolver
end

function step!(evolver::DouglasScheme, a::Vector{Float64}, t::Float64)
  t - evolver.dt > -1e-8 || error("a step towards negative time given")

  set_time!(evolver.map, max(0.0, t - evolver.dt), t)
  set_time!(evolver.bcSet, max(0.0, t - evolver.dt))

  apply_before_applying!(evolver.bcSet, evolver.map)
  y = a + evolver.dt * apply(evolver.map, a)
  apply_after_applying!(evolver.bcSet, y)

  for i = 1:get_size(evolver.map)
    rhs = y - evolver.theta * evolver.dt * apply_direction(evolver.map, i, a)
    y = solve_splitting(evolver.map, i, rhs, -evolver.theta * evolver.dt)
  end

  apply_after_solving!(evolver.bcSet, y)

  a[:] = y

  return a, evolver
end
