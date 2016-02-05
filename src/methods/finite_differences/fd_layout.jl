type FdmLinearOpLayout{I <: Integer}
  size::I
  dim::Vector{I}
  spacing::Vector{I}
end

function FdmLinearOpLayout{I <: Integer}(dim::Vector{I})
  spacing = ones(Int, length(dim))
  spacing[2:end] = cumprod(dim[1:end-1])
  sz = spacing[end] * dim[end]

  return FdmLinearOpLayout(sz, dim, spacing)
end

function get_layout_from_meshers{F1D <: Fdm1DMesher}(mesherVec::Vector{F1D})
  dim = zeros(Int, length(mesherVec))
  for i = 1:length(dim)
    dim[i] = mesherVec[i].size
  end

  return FdmLinearOpLayout(dim)
end

function neighborhood{I <: Integer}(mesherLayout::FdmLinearOpLayout, idx::I, coords::Vector{I}, i::I, offset::I)
  myIndex = idx - (coords[i] - 1) * mesherLayout.spacing[i]
  coorOffset = (coords[i] - 1) + offset

  if coorOffset < 0
    coorOffset = -coorOffset
  elseif coorOffset >= mesherLayout.dim[i]
    coorOffset = 2 * (mesherLayout.dim[i] - 1) - coorOffset
  end

  return myIndex + coorOffset * mesherLayout.spacing[i]
end

function neighborhood{I <: Integer}(mesherLayout::FdmLinearOpLayout, idx::I, coords::Vector{I}, i1::I, offset1::I, i2::I, offset2::I)
  myIndex = idx - (coords[i1] - 1) * mesherLayout.spacing[i1] - (coords[i2] - 1) * mesherLayout.spacing[i2]
  coorOffset1 = (coords[i1] - 1) + offset1

  if coorOffset1 < 0
    coorOffset1 = -coorOffset1
  elseif coorOffset1 >= mesherLayout.dim[i1]
    coorOffset1 = 2 * (mesherLayout.dim[i1] - 1) - coorOffset1
  end

  coorOffset2 = (coords[i2] - 1) + offset2
  if coorOffset2 < 0
    coorOffset2 = -coorOffset2
  elseif coorOffset2 >= mesherLayout.dim[i2]
    coorOffset2 = 2 * (mesherLayout.dim[i2] - 1) - coorOffset2
  end

  return myIndex + coorOffset1 * mesherLayout.spacing[i1] + coorOffset2 * mesherLayout.spacing[i2]
end

type FdmBoundaryConditionSet
  conditions::Vector{BoundaryCondition}
end

FdmBoundaryConditionSet() = FdmBoundaryConditionSet(Vector{BoundaryCondition}(0))

function set_time!(bcSet::FdmBoundaryConditionSet, t::Float64)
  for cond in bcSet.conditions
    set_time!(cond, t)
  end

  return bcSet
end

function apply_before_applying!(bcSet::FdmBoundaryConditionSet, op::FdmLinearOpComposite)
  for cond in bcSet.conditions
    apply_before_applying!(cond, t)
  end

  return bcSet
end

function apply_after_applying!(bcSet::FdmBoundaryConditionSet, a::Vector{Float64})
  for cond in bcSet.conditions
    apply_after_applying!(cond, t)
  end

  return bcSet
end

function apply_before_solving!(bcSet::FdmBoundaryConditionSet, op::FdmLinearOpComposite, a::Vector{Float64})
  for cond in bcSet.conditions
    apply_before_solving!(cond, op, a)
  end

  return bcSet
end

function apply_after_solving!(bcSet::FdmBoundaryConditionSet, a::Vector{Float64})
  for cond in bcSet.conditions
    apply_after_solving!(cond, t)
  end

  return bcSet
end
