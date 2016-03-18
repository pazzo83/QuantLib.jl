abstract SalvagingAlgo
type NoneSalvagingAlgo <: SalvagingAlgo end

function normalize_pseudo_root!(matrix::Matrix, pseudo::Matrix)
  sz = size(matrix)[1]
  sz == size(pseudo)[1] || error("matrix / pseudo mismatch")

  pseudoCols = size(pseudo)[2]
  # row normalization
  for i = 1:sz
    norm = 0.0
    for j = 1:pseudoCols
      norm += pseudo[i,j] * pseudo[i, j]
    end

    if norm > 0.0
      normAdj = sqrt(matrix[i, i] / norm)
      for j = 1:pseudoCols
        pseudo[i, j] *= normAdj
      end
    end
  end

  return pseudo
end

function rank_reduced_sqrt(matrix::Matrix, maxRank::Int, componentRetainedPercentage::Float64, sa::NoneSalvagingAlgo)
  sz = size(matrix)[1]
  # Symmetric Schur decomp
  jd = schurfact(matrix)
  eigenVals = jd[:values]

  # eigen vals are sorted in decreasing order
  eigenVals[end] >= -1e-16 || error("negative eigenvalue(s)")

  # factor reduction
  enough = componentRetainedPercentage * sum(eigenVals)

  if componentRetainedPercentage == 1.0
    # numerical glitches might cause some factors to be discarded
    enough *= 1.1
  end

  # retain at least one factor
  components = eigenVals[1]
  retainedFactors = 1
  i = 2
  while components < enough && i <= sz
    components += eigenVals[i]
    retainedFactors += 1
    i += 1
  end

  # output is granted to have a rank <= maxrank
  retainedFactors = min(retainedFactors, maxRank)

  diagonal = zeros(sz, retainedFactors)
  for i = 1:retainedFactors
    diagonal[i, i] = sqrt(eigenVals[i])
  end

  result = jd[:vectors] * diagonal

  normalize_pseudo_root!(matrix, result)
  
  return result
end
