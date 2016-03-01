function get_operator(process::AbstractBlackScholesProcess, grid::Vector{Float64}, residualTime::Float64, timeDependent::Bool)
  if timeDependent
    return BSMTermOperator(grid, process, residualTime)
  else
    return BSMOperator(grid, process, residualTime)
  end
end

function BSMOperator(grid::Vector{Float64}, process::AbstractBlackScholesProcess, residualTime::Float64)
  tridiag = TridiagonalOperator(length(grid))
  logGrid = LogGrid(grid)
  cc = PdeConstantCoeff(process, residualTime, state_variable(process).value)
  generate_operator!(cc, residualTime, logGrid, tridiag)

  return tridiag
end
