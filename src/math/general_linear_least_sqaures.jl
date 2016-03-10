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

  svdA = svdfact(A)
  V = svdA[:V]
  U = svdA[:U]
  w = svdA[:S]
  threshold = n * eps()
end
