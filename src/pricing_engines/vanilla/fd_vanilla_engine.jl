type FDEuropeanEngine{B <: AbstractBlackScholesProcess, I <: Integer} <: AbstractFDVanillaEngine
  process::B
  timeSteps::I
  gridPoints::I
  timeDependent::Bool
  exerciseDate::Date
  finiteDifferenceOperator::TridiagonalOperator
  intrinsicValues::SampledCurve
  BCs::Vector{BoundaryCondition}
  sMin::Float64
  center::Float64
  sMax::Float64
  prices::SampledCurve
  fdEvolverFunc::Function
end

function FDEuropeanEngine(process::AbstractBlackScholesProcess, fdEvolverFunc::Function, timeSteps::Int = 100,
                          gridPoints::Int = 100, timeDependent::Bool = false)
  prices = SampledCurve(gridPoints)
  finiteDifferenceOperator, intrinsicValues, BCs = gen_fd_vanilla_engine_params(gridPoints)
  sMin = center = sMax = 0.0
  exerciseDate = Date()

  return FDEuropeanEngine(process, timeSteps, gridPoints, timeDependent, exerciseDate, finiteDifferenceOperator,
                          intrinsicValues, BCs, sMin, center, sMax, prices, fdEvolverFunc)
end

function gen_fd_vanilla_engine_params(gridPoints::Int)
  finiteDifferenceOperator = TridiagonalOperator()
  intrinsicValues = SampledCurve(gridPoints)
  BCs = Vector{BoundaryCondition}(2)

  return finiteDifferenceOperator, intrinsicValues, BCs
end
