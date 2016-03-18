type LMMDriftCalculator
  numberOfRates::Int
  numberOfFactors::Int
  isFullFactor::Bool
  numeraire::Int
  alive::Int
  displacements::Vector{Float64}
  oneOverTaus::Vector{Float64}
  C::Matrix{Float64}
  pseudo::Matrix{Float64}
  tmp::Vector{Float64}
  e::Matrix{Float64}
  downs::Vector{Int}
  ups::Vector{Int}
end

function LMMDriftCalculator(pseudo::Matrix{Float64}, displacements::Vector{Float64}, taus::Vector{Float64}, numeraire::Int, alive::Int)
  numberOfRates = length(taus)
  prows, pcols = size(pseudo)
  numberOfFactors = pcols
  isFullFactor = numberOfFactors == numberOfRates
  tmp = zeros(numberOfRates)
  e = zeros(pcols, prows)
  downs = Vector{Int}(numberOfRates)
  ups = Vector{Int}(numberOfRates)



  # checking reqs
  numberOfRates > 0 || error("Dim out of range")
  length(displacements) == numberOfRates || error("displacements out of range")
  prows == numberOfRates || error("pseudo rows not consistent with dim")
  (pcols > 0 && pcols <= numberOfRates) || error("pseudo rows not consistent with pseudo cols")
  alive <= numberOfRates || error("Alive out of bounds")
  numeraire >= alive || error("Numeraire $(numeraire) smaller than alive $(alive)")

  # precompute 1 / taus
  oneOverTaus = 1.0 ./ taus

  # Compute covariance matrix from pseudo root
  pT = pseudo' # transpose
  C = pseudo * pT

  # Compute lower and upper extrema for (non reduced) drift calculation
  for i = alive:numberOfRates
    downs[i] = min(i+1, numeraire)
    ups[i] = max(i+1, numeraire)
  end

  return LMMDriftCalculator(numberOfRates, numberOfFactors, isFullFactor, numeraire, alive, displacements, oneOverTaus, C, pseudo,
                            tmp, e, downs, ups)
end

function compute_reduced!(lmm::LMMDriftCalculator, forwards::Vector{Float64}, drifts::Vector{Float64})
  # Compute drifts with factor reduction using pseudo square root of covariance matrix

  # precompute forwards factor
  for i = lmm.alive:lmm.numberOfRates
    lmm.tmp[i] = (forwards[i] + lmm.displacements[i]) / (lmm.oneOverTaus[i] + forwards[i])
  end

  # Enforce initialization
  for r = 1:lmm.numberOfFactors
    lmm.e[r, max(1, lmm.numeraire-1)] = 0.0
  end
  # Now compute drifts: take the numeraire P_N (numeraire = N) as the reference point
  # divide the summation into 3 steps

  # 1st Step: the drift corresponding to the numeraire P_N is zero
  # (if N=1 no drift is null, if N = numberOfRates the last drift is null)
  if lmm.numeraire > 1
    drifts[lmm.numeraire - 1] = 0.0
  end

  # 2nd Step: then, move backwards from N-2 (included) back to alive(included)
  # (if N=1, jump to 3rd Step, if N = numberOfRates, the e[r, N-1] are correctly initialized)
  for i = lmm.numeraire - 2:-1:lmm.alive
    drifts[i] = 0.0
    for r = 1:lmm.numberOfFactors
      lmm.e[r, i] = lmm.e[r, i+1] + lmm.tmp[i+1] * lmm.pseudo[i+1, r]
      drifts[i] -= lmm.e[r, i] * lmm.pseudo[i, r]
    end
  end

  # 3rd step: now move forward from N (included) up to n (excluded)
  # (if N = 0 this is the only relevant computation)
  for i = lmm.numeraire:lmm.numberOfRates
    drifts[i] = 0.0
    for r = 1:lmm.numberOfFactors
      if i == 1
        lmm.e[r, i] = lmm.tmp[i] * lmm.pseudo[i, r]
      else
        lmm.e[r, i] = lmm.e[r, i-1] + lmm.tmp[i] * lmm.pseudo[i, r]
      end
      drifts[i] += lmm.e[r, i] * lmm.pseudo[i, r]
    end
  end
  return lmm, drifts
end

function compute!(lmm::LMMDriftCalculator, forwards::Vector{Float64}, drifts::Vector{Float64})
  # TODO check for full factor
  compute_reduced!(lmm, forwards, drifts)

  return lmm, drifts
end
