type SwapForwardMappings end

function coterminal_swap_forward_jacobian(::SwapForwardMappings, cs::CurveState)
  n = cs.numberOfRates
  f = cs.forwardRates
  tau = cs.rateTaus

  a = Vector{Float64}(n)
  for k in eachindex(a)
    a[k] = discount_ratio(cs, k, n) - 1.0
  end

  jacobian = zeros(n, n)
  for i = 1:n, j = i:n
    bi = coterminal_swap_annuity(cs, n, i)
    bj = coterminal_swap_annuity(cs, n, j)
    jacobian[i, j] = tau[j] / coterminal_swap_annuity(cs, j+1, i) + tau[j] / (1.0 + f[j] * tau[j]) * (-a[j] * bi + a[i] * bj) / (bi * bi)
  end

  return jacobian
end


function coterminal_swap_zed_matrix(sfm::SwapForwardMappings, cs::CurveState, displacement::Float64)
  n = cs.numberOfRates
  zMatrix = coterminal_swap_forward_jacobian(sfm, cs)

  f = cs.forwardRates
  sr = coterminal_swap_rate(cs)

  for i = 1:n, j = i:n
    zMatrix[i, j] *= (f[j] + displacement) / (sr[i] + displacement)
  end

  return zMatrix    
end
