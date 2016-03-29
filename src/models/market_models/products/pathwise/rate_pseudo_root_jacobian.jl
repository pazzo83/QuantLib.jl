type RatePseudoRootJacobianAllElements
  pseudoRoot::Matrix{Float64}
  aliveIndex::Int
  taus::Vector{Float64}
  pseudoBumps::Vector{Matrix{Float64}}
  displacements::Vector{Float64}
  factors::Int
  e::Matrix{Float64}
  ratios::Vector{Float64}
end

function RatePseudoRootJacobianAllElements(pseudoRoot::Matrix{Float64},
                                            aliveIndex::Int,
                                            numeraire::Int,
                                            taus::Vector{Float64},
                                            displacements::Vector{Float64})
  pr_rows, factors = size(pseudoRoot)
  e = Matrix{Float64}(pr_rows, factors)
  numberRates = length(taus)
  ratios = Vector{Float64}(numberRates)

  aliveIndex == numeraire || error("we can only do discrete compounding MM account so alive Index must equal numeraire")
  pr_rows == numberRates || error("pseudoRoot rows != taus size")
  length(displacements) == numberRates || error("displacements size != taus size")

  return RatePseudoRootJacobianAllElements(pseudoRoot, aliveIndex, taus, Vector{Matrix{Float64}}(), displacements, factors, e, ratios)
end

function get_bumps!(ratePseudo::RatePseudoRootJacobianAllElements, oldRates::Vector{Float64}, discountRatios::Vector{Float64}, newRates::Vector{Float64},
                    gaussians::Vector{Float64}, B::Vector{Matrix{Float64}})
  numberRates = length(ratePseudo.taus)

  length(B) == numberRates || error("B must be the same size as numberRates")
  for j in eachindex(B)
    size(B) == (numberRates, ratePseudo.factors) || error("matrix in B not of proper size")
  end

  for j = ratePseudo.aliveIndex:numberRates
    ratePseudo.ratios[j] = (oldRates[j] + ratePseudo.displacements[j]) * discountRatios[j+1]
  end

  for f = 1:factors
    ratePseudo.e[ratePseudo.aliveIndex, f] = 0.0
    for j = ratePseudo.aliveIndex+1:numberRates
      ratePseudo.e[j, f] = ratePseudo.e[j-1, f] + ratePseudo.ratios[j-1] * ratePseudo.pseudoBumps[j-1, f]
    end
  end

  # nullify B for rates that have already reset
  for j = 1:ratePseudo.aliveIndex, k = 1:numberRates, f = 1:ratePseudo.factors
    B[j][k, f] = 0.0
  end

  for f = 1:ratePseudo.factors, j = ratePseudo.aliveIndex:numberRates
    for k = ratePseudo.aliveIndex:j
      B[j][k, f] = newRates[j] * ratePseudo.ratios[k] * ratePseudo.taus[k] * ratePseudo.pseudoRoot[j, f]
    end

    tmp = 2.0 * ratePseudo.ratios[j] * ratePseudo.taus[j] * ratePseudo.pseudoRoot[j, f]
    tmp -= ratePseudo.pseudoRoot[j, f]
    tmp += ratePseudo.e[j, f] * ratePseudo.taus[j]
    tmp += gaussians[f]
    tmp *= (newRates[j] + ratePseudo.displacements[j])

    B[j][j, f] = tmp

    for k = 1:ratePseudo.aliveIndex
      B[j][k, f] = 0.0
    end

    for k = j+1:numberRates
      B[j][k, f] = 0.0
    end
  end

  return ratePseudo, B
end
