# Solvers
abstract Solver1D

type SolverInfo
  maxEvals::Int
  lowerBoundEnforced::Bool
  upperBoundEnforced::Bool
  lowerBound::Float64
  upperBound::Float64
end

# solver functions
function solve(solver::Solver1D, f::Function, accuracy::Float64, guess::Float64, step::Float64)
  ## This method returns the 0 of a function determined by a given accuracy.  This method using bracketing
  ## routine to which an intial guess must be supplied as well as a step
  growth_factor = 1.6
  flipflop = -1

  root = guess
  fxMax = f(root)

  # monotonically crescent bias, as in optionValue(volatility)
  if is_close(fxMax, 0.0)
    return root
  elseif fxMax > 0.0
    xMin = enforced_bounds(solver, root - step)
    fxMin = f(xMin)
    xMax = root
  else
    xMin = root
    fxMin = fxMax
    xMax = enforced_bounds(solver, root + step)
    fxMax = f(xMax)
  end

  eval_num = 2

  while eval_num < solver.solverInfo.maxEvals
    if fxMin * fxMax <= 0.0
      if is_close(fxMin, 0.0)
        return xMin
      end
      if is_close(fxMax, 0.0)
        return xMax
      end
      root = (xMax + xMin) / 2.0
      return _solve(solver, f, accuracy, xMin, xMax, fxMin, fxMax, root, eval_num)
    end

    if abs(fxMin) < abs(fxMax)
      xMin = enforced_bounds(solver, xMin + growth_factor * (xMin - xMax))
      fxMin = f(xMin)
    elseif abs(fxMin) > abs(fxMax)
      xMax = enforced_bounds(solver, xMax + growth_factor * (xMax - xMin))
      fxMax = f(xMax)
    elseif flipflop == -1
      xMin = enforced_bounds(solver, xMin + growth_factor * (xMin - xMax))
      fxMin = f(xMin)
      eval_num += 1
      flipflop = 1
    elseif flipflop == 1
      xMax = enforced_bounds(solver, xMax + growth_factor * (xMax - xMin))
      fxMax = f(xMax)
      flipflop = -1
    end

    eval_num += 1
  end

  return root
end

function solve(solver::Solver1D, f::Function, accuracy::Float64, guess::Float64, xMin::Float64, xMax::Float64)
  fxMin = f(xMin)
  if is_close(fxMin, 0.0)
    return xMin
  end

  fxMax = f(xMax)
  if is_close(fxMax, 0.0)
    return xMax
  end

  eval_num = 2

  return _solve(solver, f, accuracy, xMin, xMax, fxMin, fxMax, guess, eval_num)
end

function enforced_bounds{S <: Solver1D}(solver::S, x::Float64)
  if solver.solverInfo.lowerBoundEnforced && x < solver.solverInfo.lowerBound
    return solver.solverInfo.lowerBound
  end

  if solver.solverInfo.upperBoundEnforced && x > solver.solverInfo.upperBound
    return solver.solverInfo.upperBound
  end

  return x
end
