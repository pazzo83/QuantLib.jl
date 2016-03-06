# Thank you Distributions.jl!
using Distributions
# normal distribution methods
function distribution_derivative(w::Normal, x::Float64)
  _sigma = std(w)
  xn = (x - mean(w)) / _sigma
  return pdf(w, xn) / _sigma
end

function peizer_pratt_method_2_inversion(z::Float64, n::Int)
  isodd(n) || error("n must be an odd number, $n not allowed")

  res = (z / (n + 1.0 / 3.0 + 0.1 / (n + 1.0)))
  res *= res
  res = exp(-res * (n + 1.0 / 6.0))
  res = 0.5 + (z > 0.0 ? 1.0 : -1.0) * sqrt(0.25 * (1.0 - res))
  return res
end
