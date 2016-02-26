@enum Side NoneSide UpperSide LowerSide

type DefaultBoundaryCondition <: BoundaryConditionType end
type NeumannBC <: BoundaryConditionType end
type DirichletBC <: BoundaryConditionType end

type BoundaryCondition{B <: BoundaryConditionType}
  condType::B
  value::Float64
  side::Side
end
