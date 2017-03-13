abstract SalvagingAlgo
type NoneSalvagingAlgo <: SalvagingAlgo end

function normalize_pseudo_root!(matrix::Matrix, pseudo::Matrix)
  sz = size(matrix)[1]
  sz == size(pseudo)[1] || error("matrix / pseudo mismatch")

  pseudoCols = size(pseudo)[2]
  # row normalization
  @inbounds @simd for i = 1:sz
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
    @inbounds components += eigenVals[i]
    retainedFactors += 1
    i += 1
  end

  # output is granted to have a rank <= maxrank
  retainedFactors = min(retainedFactors, maxRank)

  diagonal = zeros(sz, retainedFactors)
  @simd for i = 1:retainedFactors
    @inbounds diagonal[i, i] = sqrt(eigenVals[i])
  end

  result = jd[:vectors] * diagonal

  normalize_pseudo_root!(matrix, result)

  return result
end

function ql_norm_squared(v::Matrix{Float64}, row::Int)
  x = 0.0
  @simd for i = 1:size(v)[2]
    @inbounds x += v[row, i] * v[row, i]
  end

  return x
end

function ql_norm(v::Matrix{Float64}, row::Int)
  return sqrt(ql_norm_squared(v, row))
end

function inner_product(v::Matrix{Float64}, row1::Int, w::Matrix{Float64}, row2::Int)
  x = 0.0
  @simd for i = 1:size(v)[2]
    @inbounds x += v[row1, i] * w[row2, i]
  end

  return x
end

type OrthogonalProjection
  originalVectors::Matrix{Float64}
  multiplierCutoff::Float64
  numberVectors::Int
  numberValidVectors::Int
  dimension::Int
  validVectors::BitArray{1}
  projectedVectors::Vector{Vector{Float64}}
  orthoNormalizedVectors::Matrix{Float64}
end

function OrthogonalProjection(originalVectors::Matrix{Float64}, multiplierCutoff::Float64, tolerance::Float64)
  numberVectors, dimension = size(originalVectors)
  validVectors = trues(numberVectors)
  orthoNormalizedVectors = Matrix{Float64}(numberVectors, dimension)
  projectedVectors = Vector{Vector{Float64}}()

  currentVector = Vector{Float64}(dimension)
  @inbounds @simd for j in eachindex(validVectors)
    if validVectors[j]
      for k = 1:numberVectors
        for m = 1:dimension
          orthoNormalizedVectors[k, m] = originalVectors[k, m]
        end

        if k != j && validVectors[k]
          for l = 1:k
            if validVectors[l] && l != j
              dotProduct = inner_product(orthoNormalizedVectors, k, orthoNormalizedVectors, l)
              for n = 1:dimension
                orthoNormalizedVectors[k, n] -= dotProduct * orthoNormalizedVectors[l, n]
              end
            end
          end

          normBeforeScaling = ql_norm(orthoNormalizedVectors, k)

          if normBeforeScaling < tolerance
            validVectors[k] = false
          else
            normBeforeScalingRecip = 1.0/normBeforeScaling
            for m = 1:dimension
              orthoNormalizedVectors[k, m] *= normBeforeScalingRecip
            end
          end # end of else (norm < tolerance)
        end # end of if k != j && validVectors[k]
      end # end of for k = 1:numberVectors
      # we now have an o.n. basis for everything except j
      prevNormSquared = ql_norm_squared(originalVectors, j)

      for r = 1:numberVectors
        if validVectors[r] && r != j
          dotProduct = inner_product(orthoNormalizedVectors, j, orthoNormalizedVectors, r)
          for s = 1:dimension
            orthoNormalizedVectors[j, s] -= dotProduct * orthoNormalizedVectors[r, s]
          end
        end
      end

      projectionOnOriginalDirection = inner_product(originalVectors, j, orthoNormalizedVectors, j)
      sizeMultiplier = prevNormSquared / projectionOnOriginalDirection

      if abs(sizeMultiplier) < multiplierCutoff
        for t = 1:dimension
          currentVector[t] = orthoNormalizedVectors[j, t] * sizeMultiplier
        end
      else
        validVectors[j] = false
      end
    end # end of validVectors[j]

    push!(projectedVectors, copy(currentVector))
  end # end of j loop

  numberValidVectors = 0
  @simd for i in eachindex(validVectors)
    @inbounds numberValidVectors += validVectors[i] ? 1 : 0
  end

  return OrthogonalProjection(originalVectors, multiplierCutoff, numberVectors, numberValidVectors, dimension, validVectors, projectedVectors, orthoNormalizedVectors)
end

get_vector(op::OrthogonalProjection, i::Int) = op.projectedVectors[i]
