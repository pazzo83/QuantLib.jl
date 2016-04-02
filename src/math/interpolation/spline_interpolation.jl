using Dierckx

type BicubicSpline
  spline::Dierckx.Spline2D
end

BicubicSpline{T <: Real}(x::Vector{T}, y::Vector{T}, z::Matrix{T}) = BicubicSpline(Dierckx.Spline2D(x, y, z))

type NaturalCubicSpline{T <: Number} <: Interpolation
  x_vert::Vector{T}
  y_vert::Vector{T}
  b::Vector{T}
  c::Vector{T}
  d::Vector{T}
end

# adapted from here: http://sepwww.stanford.edu/sep/sergey/128A/answers6.pdf
function NaturalCubicSpline{T <: Real}(x_vert::Vector{T}, y_vert::Vector{T})
  n = length(x_vert)
  h = zeros(n - 1)
  b = zeros(n - 1)
  a = zeros(n - 2)
  g = zeros(n - 2)
  c = zeros(n)
  d = zeros(n - 1)

  for i = 1:n - 1
    h[i] = x_vert[i + 1] - x_vert[i]
    b[i] = (y_vert[i + 1] - y_vert[i]) / h[i]
  end

  for i = 1:n - 2
    a[i] = 2.0 * (h[i] + h[i + 1])
    # u[i] = 6 * (b[i + 1] - b[i])
    g[i] = b[i + 1] - b[i]
  end

  Alu = lufact(Tridiagonal(h[2:end-1], a[1:end], h[2:end-1]))
  c[2:end - 1] = Alu \ g

  for i=1:n - 1
    d[i] = (c[i+1] - c[i])/h[i]
    b[i] -= (2.0 * c[i] + c[i + 1]) * h[i]
    c[i] *= 3.0
  end

  return NaturalCubicSpline(x_vert, y_vert, b, c, d)
end

function call{T <: Real}(spl::NaturalCubicSpline, my_x::T)
  # get coefficients
  # b, c, d = gen_splines(x_vert, y_vert)

  x_idx = searchsortedlast(spl.x_vert, my_x)

  if x_idx == 0
    x_idx = 1
  elseif x_idx == length(spl.x_vert)
    x_idx = length(spl.x_vert) - 1
  end


  diff = my_x - spl.x_vert[x_idx]

  return spl.y_vert[x_idx] + diff*(spl.b[x_idx] + diff * (spl.c[x_idx] + (diff * spl.d[x_idx])))
end

## Cubic Interpolation ##
type CubicInterpolation{D <: DerivativeApprox, B1 <: BoundaryCondition, B2 <: BoundaryCondition} <: Interpolation
  derivativeApprox::D
  leftBoundaryCondition::B1
  rightBoundaryCondition::B2
  leftValue::Float64
  rightValue::Float64
  monotonic::Bool
  x_vals::Vector{Float64}
  y_vals::Vector{Float64}
  a::Vector{Float64}
  b::Vector{Float64}
  c::Vector{Float64}
  tmp::Vector{Float64}
  dx::Vector{Float64}
  S::Vector{Float64}
  n::Int
  L::TridiagonalOperator
end

function CubicInterpolation(dApprox::DerivativeApprox, leftBoundary::BoundaryCondition, rightBoundary::BoundaryCondition, x_vals::Vector{Float64}, y_vals::Vector{Float64},
                            leftValue::Float64 = 0.0, rightValue::Float64 = 0.0, monotonic::Bool = true)
  n = length(x_vals)
  ci = CubicInterpolation(dApprox, leftBoundary, rightBoundary, leftValue, rightValue, monotonic, x_vals, y_vals, zeros(n - 1), zeros(n - 1), zeros(n - 1), zeros(n), zeros(n - 1), zeros(n - 1),
                          n, TridiagonalOperator(n))

  update!(ci)
  return ci
end

call(interp::CubicInterpolation, x::Float64) = value(interp, x)

# type aliases #
typealias SplineCubicInterpolation{D, B1, B2} CubicInterpolation{Spline, B1, B2} # First derivative approximation

function cubicInterpolationPolynomialDerivative(x_vals::Vector{Float64}, y_vals::Vector{Float64}, x::Float64)
  a, b, c, d = x_vals
  u, v, w, z = y_vals

  return (-((((a-c)*(b-c)*(c-x)*z-(a-d)*(b-d)*(d-x)*w)*(a-x+b-x)+((a-c)*(b-c)*z-(a-d)*(b-d)*w)*(a-x)*(b-x))*(a-b)+
          ((a-c)*(a-d)*v-(b-c)*(b-d)*u)*(c-d)*(c-x)*(d-x)+((a-c)*(a-d)*(a-x)*v-(b-c)*(b-d)*(b-x)*u)
          *(c-x+d-x)*(c-d)))/((a-b)*(a-c)*(a-d)*(b-c)*(b-d)*(c-d))
end

function update_left_boundary_condition!(::Lagrange, interp::SplineCubicInterpolation)
  set_first_row!(interp.L, 1.0, 0.0)
  interp.tmp[1] = cubicInterpolationPolynomialDerivative(interp.x_vals[1:4], interp.y_vals[1:4], interp.x_vals[1])

  return interp
end

function update_right_boundary_condition!(::Lagrange, interp::SplineCubicInterpolation)
  set_last_row!(interp.L, 0.0, 1.0)
  interp.tmp[interp.n] = cubicInterpolationPolynomialDerivative(interp.x_vals[interp.n - 3:interp.n], interp.y_vals[interp.n - 3:interp.n], interp.x_vals[interp.n])

  return interp
end

function update!(interp::SplineCubicInterpolation)
  n = interp.n

  for i = 1:n - 1
    interp.dx[i] = interp.x_vals[i + 1] - interp.x_vals[i]
    interp.S[i] = (interp.y_vals[i + 1] - interp.y_vals[i]) / interp.dx[i]
  end

  for i = 2:n - 1
    set_mid_row!(interp.L, i, interp.dx[i], 2.0 * (interp.dx[i] + interp.dx[i - 1]), interp.dx[i - 1])
    interp.tmp[i] = 3.0 * (interp.dx[i] * interp.S[i - 1] + interp.dx[i - 1] * interp.S[i])
  end

  # update left boundary
  update_left_boundary_condition!(interp.leftBoundaryCondition, interp)

  # update right boundary
  update_right_boundary_condition!(interp.rightBoundaryCondition, interp)

  # solve
  solve_for!(interp.L, interp.tmp, interp.tmp)
  # interp.tmp = solve_for(interp.L, interp.tmp)

  if interp.monotonic
    corretion = 0.0
    pm = pu = pd = M = 0.0
    for i = 1:n
      if i == 1
        if interp.tmp[i] * interp.S[i] > 0.0
          correction = interp.tmp[i] / abs(interp.tmp[i]) * min(abs(interp.tmp[i]), abs(3.0 * interp.S[1]))
        else
          correction = 0.0
        end

        if correction != interp.tmp[i]
          interp.tmp[i] = correction
        end
      elseif i == n
        if interp.tmp[i] * interp.S[n - 1] > 0.0
          correction = interp.tmp[i] / abs(interp.tmp[i]) * min(abs(interp.tmp[i]), abs(3.0 * interp.S[n - 1]))
        else
          correction = 0.0
        end

        if correction != interp.tmp[i]
          interp.tmp[i] = correction
        end
      else
        pm = (interp.S[i - 1] * interp.dx[i] + interp.S[i] * interp.dx[i - 1]) / (interp.dx[i - 1] + interp.dx[i])
        M = 3.0 * min(min(abs(interp.S[i - 1]), abs(interp.S[i])), abs(pm))

        if i > 2
          if (interp.S[i - 1] - interp.S[i - 2]) * (interp.S[i] - interp.S[i - 1]) > 0.0
            pd = (interp.S[i - 1] * (2.0 * interp.dx[i - 1] + interp.dx[i - 2]) - interp.S[i - 2] * interp.dx[i - 1]) / (interp.dx[i - 2] + interp.dx[i - 1])

            if pm * pd > 0.0 && pm * (interp.S[i - 1] - interp.S[i - 2]) > 0.0
              M = max(M, 15 * min(abs(pm), abs(pd)))
            end
          end
        end
        if i < n - 1
          if (interp.S[i] - interp.S[i - 1]) * (interp.S[i + 1] - interp.S[i]) > 0.0
            pu = (interp.S[i] * (2.0 * interp.dx[i] + interp.dx[i + 1]) - interp.S[i + 1] * interp.dx[i]) / (interp.dx[i] + interp.dx[i + 1])
            if pm * pu > 0.0 && -pm * (interp.S[i] - interp.S[i - 1]) > 0.0
              M = max(M, 1.5 * min(abs(pm), abs(pu)))
            end
          end
        end

        if interp.tmp[i] * pm > 0.0
          correction = interp.tmp[i] / abs(interp.tmp[i]) * min(abs(interp.tmp[i]), M)
        else
          correction = 0.0
        end

        if correction != interp.tmp[i]
          interp.tmp[i] = correction
        end
      end
    end
  end

  # coefficients
  for i = 1:interp.n - 1
    interp.a[i] = interp.tmp[i]
    interp.b[i] = (3.0 * interp.S[i] - interp.tmp[i + 1] - 2.0 * interp.tmp[i]) / interp.dx[i]
    interp.c[i] = (interp.tmp[i + 1] + interp.tmp[i] - 2.0 * interp.S[i]) / (interp.dx[i] * interp.dx[i])
  end

  return interp
end

function value(interp::CubicInterpolation, x::Float64)
  j = locate(interp, x)
  dx = x - interp.x_vals[j]
  return interp.y_vals[j] + dx * (interp.a[j] + dx * (interp.b[j] + dx * interp.c[j]))
end
