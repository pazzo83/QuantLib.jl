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
