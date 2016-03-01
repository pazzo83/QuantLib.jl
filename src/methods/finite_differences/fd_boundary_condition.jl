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

NeumannBC(value::Float64, side::BC_Side) = BoundaryCondition(NeumannBCType(), value, side)
