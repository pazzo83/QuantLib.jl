type NewtonSolver <: Solver1D
  solverInfo::SolverInfo
end
NewtonSolver{I <: Integer}(maxEvals::I = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0) =
  NewtonSolver(SolverInfo(maxEvals, lowerBoundEnforced, upperBoundEnforced, lowerBound, upperBound))

# Newton Solver (safe)
function _solve{I <: Integer}(solver::NewtonSolver, f::Function, accuracy::Float64, xMin::Float64, xMax::Float64, fxMin::Float64, fxMax::Float64, root::Float64, eval_num::I)
  # Orient the search so that f(xl) < 0
  if fxMin < 0.0
    xl = xMin
    xh = xMax
  else
    xh = xMin
    xl = xMax
  end

  # the stepsize before last
  dxold = xMax - xMin

  # and the last step
  dx = dxold

  froot = f(root)
  dfroot = f(root, Derivative()) # derivative, hurray multiple dispatch!
  eval_num += 1

  while eval_num < solver.solverInfo.maxEvals
    # bisect if out of range or not decreasing fast enough
    if ((root - xh) * dfroot - froot) * ((root - xl) * dfroot - froot) > 0.0 || abs(2.0 * froot) > abs(dxold * dfroot)
      dxold = dx
      dx = (xh - xl) / 2.0
      root = xl + dx
    else
      dxold = dx
      dx = froot / dfroot
      root -= dx
    end

    # Convergence criterion
    if abs(dx) < accuracy
      # root = f(root) # check this
      eval_num += 1
      return root
    end

    froot = f(root)
    dfroot = f(root, Derivative())
    eval_num += 1
    if froot < 0.0
      xl = root
    else
      xh = root
    end
  end

  error("Maximum number of function evals exceeded!")
end
