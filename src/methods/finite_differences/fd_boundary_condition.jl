# @enum Side NoneSide UpperSide LowerSide
type NoneSide <: BC_Side end
type UpperSide <: BC_Side end
type LowerSide <: BC_Side end

type DefaultBoundaryCondition <: BoundaryConditionType end
type NeumannBCType <: BoundaryConditionType end
type DirichletBCType <: BoundaryConditionType end

type BoundaryCondition{B <: BoundaryConditionType, S <: BC_Side}
  condType::B
  value::Float64
  side::S
end

typealias NeumannBC{S} BoundaryCondition{NeumannBCType, S}
build_NeumannBC(value::Float64, side::BC_Side) = BoundaryCondition(NeumannBCType(), value, side)

set_time!(bc::NeumannBC, ::Float64) = bc # returns self

apply_before_applying!(bc::NeumannBC, L::TridiagonalOperator) = apply_before_applying!(bc, bc.side, L)

function apply_before_applying!(bc::NeumannBC, ::LowerSide, L::TridiagonalOperator)
  set_first_row!(L, -1.0, 1.0)
  return L
end

function apply_before_applying!(bc::NeumannBC, ::UpperSide, L::TridiagonalOperator)
  set_last_row!(L, -1.0, 1.0)
  return L
end

apply_after_applying!(bc::NeumannBC, u::Vector{Float64}) = apply_after_applying!(bc, bc.side, u)

function apply_after_applying!(bc::NeumannBC, ::LowerSide, u::Vector{Float64})
  u[1] = u[2] - bc.value
  return u
end

function apply_after_applying!(bc::NeumannBC, ::UpperSide, u::Vector{Float64})
  u[end] = u[end-1] + bc.value
  return u
end

apply_before_solving!(bc::NeumannBC, L::TridiagonalOperator, rhs::Vector{Float64}) = apply_before_solving!(bc, bc.side, L, rhs)

function apply_before_solving!(bc::NeumannBC, ::LowerSide, L::TridiagonalOperator, rhs::Vector{Float64})
  set_first_row!(L, -1.0, 1.0)
  rhs[1] = bc.value
  return L, rhs
end

function apply_before_solving!(bc::NeumannBC, ::UpperSide, L::TridiagonalOperator, rhs::Vector{Float64})
  set_last_row!(L, -1.0, 1.0)
  rhs[end] = bc.value
  return L, rhs
end

apply_after_solving!(bc::NeumannBC, ::Vector{Float64}) = bc # returns self
