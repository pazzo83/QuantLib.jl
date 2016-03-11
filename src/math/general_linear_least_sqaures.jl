type GeneralLinearLeastSquares
  a::Vector{Float64}
  err::Vector{Float64}
  residuals::Vector{Float64}
  standardErrors::Vector{Float64}
end

function GeneralLinearLeastSquares{T}(x::Vector{Float64}, y::Vector{Float64}, v::Vector{T})
  n = length(v)

  a = zeros(n)
  err = zeros(n)
  residuals = zeros(length(y))
  standardErrors = zeros(n)
  glls = GeneralLinearLeastSquares(a, err, residuals, standardErrors)
  calculate!(glls, x, y, v)

  return glls
end

function calculate!{T}(glls::GeneralLinearLeastSquares, x::Vector{Float64}, y::Vector{Float64}, v::Vector{T})
  n = length(glls.residuals)
  m = length(glls.err)

  A = zeros(n, m)

  for i in eachindex(v)
    A[:, i] = map(v[i], x)
  end
  # println(A[:, 2])
  # println("-----------------")
  # println(A[:, 3])
  # println("-----------------")
  # println(A[:, 4])
  # error("DIE")
  svdA = svdfact(A)
  V = svdA[:Vt] * -1.0
  U = svdA[:U] * -1.0
  w = svdA[:S]
  # svdA = SVD(A)
  # V = svdA.V
  # U = svdA.U
  # w = svdA.s
  threshold = n * eps()
  for i in eachindex(w)
    if w[i] > threshold
      u = dot(U[:, i], y) / w[i]
      for j in eachindex(glls.a)
        glls.a[j] += u * V[j, i]
        glls.err[j] += V[j, i] * V[j, i] / (w[i] * w[i])
      end
    end
  end
  # println(glls.a)
  # error("DIE")
  glls.err = sqrt(glls.err)
  tmp = A * glls.a
  glls.residuals = tmp - y

  chiSq = dot(glls.residuals, glls.residuals)

  glls.standardErrors = glls.err * sqrt(chiSq / (n - 2))

  return glls
end

get_coefficients(glls::GeneralLinearLeastSquares) = glls.a
