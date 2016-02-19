type FiniteDifferenceNewtonSafe <: Solver1D
  solverInfo::SolverInfo
end
FiniteDifferenceNewtonSafe{I <: Integer}(maxEvals::I = 100, lowerBoundEnforced::Bool = false, upperBoundEnforced::Bool = false, lowerBound::Float64 = 0.0, upperBound::Float64 = 0.0) =
  FiniteDifferenceNewtonSafe(SolverInfo(maxEvals, lowerBoundEnforced, upperBoundEnforced, lowerBound, upperBound))

# Finite Differences Solver
function _solve{I <: Integer}(solver::FiniteDifferenceNewtonSafe, f::Function, accuracy::Float64, xMin::Float64, xMax::Float64, fxMin::Float64, fxMax::Float64,
                root::Float64, eval_num::I)

  max_evals = solver.solverInfo.maxEvals
  # orienting search such that f(xl) < 0
  if fxMin < 0.0
    xl = xMin
    xh = xMax
  else
    xh = xMin
    xl = xMax
  end

  froot = f(root)
  eval_num += 1

  # first order finite difference derivative
  dfroot = xMax - root < root - xMin ? (fxMax - froot) / (xMax - root) : (fxMin - froot) / (xMin - root)

  dx = xMax - xMin
  while eval_num <= max_evals
    frootold = froot
    rootold = root
    dxold = dx
    # bisect if out of range or not decreasing fast enough
    if (((root - xh) * dfroot - froot) * ((root - xl) * dfroot - froot)) > 0.0 || abs(2.0 * froot) > abs(dxold * dfroot)
      dx = (xh - xl) / 2.0
      root = xl + dx
      # if root estimate computed is close to previous one, calculate dfroot at root ane xh rather than
      # root and rootold
      if is_close(root, rootold, 2500)
        rootold = xh
        frootold = f(xh)
      end
    else
      # Newton
      dx = froot / dfroot
      root -= dx
    end

    # Convergence check
    if abs(dx) < accuracy
      return root
    end

    froot = f(root)
    eval_num += 1
    dfroot = (frootold - froot) / (rootold - root)

    if (froot < 0.0)
      xl = root
    else
      xh = root
    end
  end

  error("Maximum number of function evals exceeded!")
end
