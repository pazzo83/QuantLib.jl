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

type MixedScheme <: FdScheme
  L::TridiagonalOperator
  Ident::TridiagonalOperator
  dt::Float64
  theta::Float64
  bcSet::FdmBoundaryConditionSet
  explicitPart::TridiagonalOperator
  implicitPart::TridiagonalOperator

  function MixedScheme(L::TridiagonalOperator, theta::Float64, bcSet::FdmBoundaryConditionSet)
    Ident = TridiagIdentity(L.n)
    dt = 0.0
    return new(L, Ident, dt, theta, bcSet)
  end
end

# add'l constructors
CrankNelson(L::TridiagonalOperator, bcs::FdmBoundaryConditionSet) = MixedScheme(L, 0.5, bcs)

function set_step!(evolver::MixedScheme, dt::Float64)
  evolver.dt = dt
  if evolver.theta != 1.0 # there is an explicit part
    evolver.explicitPart = evolver.Ident - ((1.0 - evolver.theta) * evolver.dt) * evolver.L
  end

  if evolver.theta != 0.0 # there is an implicit part
    evolver.implicitPart = evolver.Ident + (evolver.theta * evolver.dt) * evolver.L
  end

  return evolver
end

function step!(evolver::MixedScheme, a::Vector{Float64}, t::Float64)
  for i in eachindex(evolver.bcSet.conditions)
    set_time!(evolver.bcSet.conditions[i], t)
  end

  if evolver.theta != 1.0 # there is an explicit part
    # TODO Tridiagonal time dependent
    for i in eachindex(evolver.bcSet.conditions)
      apply_before_applying!(evolver.bcSet.conditions[i], evolver.explicitPart)
    end
    a[:] = apply_to(evolver.explicitPart, a)
    for i in eachindex(evolver.bcSet.conditions)
      apply_after_applying!(evolver.bcSet.conditions[i], a)
    end
  end

  if evolver.theta != 0.0 # there is an implicit part
    # TODO Tridiagonal time dependent
    for i in eachindex(evolver.bcSet.conditions)
      apply_before_solving!(evolver.bcSet.conditions[i], evolver.implicitPart, a)
    end
    solve_for!(evolver.implicitPart, a, a)
    for i in eachindex(evolver.bcSet.conditions)
      apply_after_solving!(evolver.bcSet.conditions[i], a)
    end
  end

  return a, evolver
end

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
