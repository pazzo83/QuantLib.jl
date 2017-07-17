struct TimeHomogeneousForwardCorrelation end

function evolved_matrices(::Type{TimeHomogeneousForwardCorrelation}, fwdCorrelation::Matrix{Float64})
  numberOfRates = size(fwdCorrelation)[1]
  correlations = Matrix{Float64}[zeros(numberOfRates, numberOfRates) for i = 1:numberOfRates]
  # correlations = fill(zeros(numberOfRates, numberOfRates), numberOfRates)

  @inbounds @simd for k in eachindex(correlations)
    # proper diagonal values
    for i = k:numberOfRates
      correlations[k][i, i] = 1.0
    end
    # copy only time homogeneous values
    for i = k:numberOfRates, j = k:i
      correlations[k][i,j] = correlations[k][j,i] = fwdCorrelation[i-(k-1), j-(k-1)]
    end
  end
  return correlations
end
