function generic_longstaff_schwartz_regression!(simulationData::Vector{Vector{NodeData}}, basisCoefficients::Vector{Vector{Float64}})
  steps = length(simulationData)
  resize!(basisCoefficients, steps-1)

  for i = steps:-1:2
    exerciseData = simulationData[i]

    # 1) Find the covariance matrix of basis function values and deflated cash-flows
    N = length(exerciseData[1].values)
    temp = Vector{Float64}(N+1)
    stats = GenericSequenceStats(N+1, length(exerciseData))

    for j in eachindex(exerciseData)
      if exerciseData[j].isValid
        temp[1:end-1] = exerciseData[j].values
        temp[end] = exerciseData[j].cumulatedCashFlows - exerciseData[j].controlValue

        add_sample!(stats, temp, j)
      end
    end

    means = stats_mean(stats)
    covMat = stats_covariance(stats)

    C = Matrix{Float64}(N, N)
    target = Vector{Float64}(N)

    for k in eachindex(target)
      target[k] = covMat[k, N+1] + means[k] * means[N+1]
      for l = 1:k
        C[k, l] = C[l, k] = covMat[k,l] + means[k] * means[l]
      end
    end

    # 2) Solve for least squares regression
    alphas = svdfact(C) \ target
    # resize!(basisCoefficients[i-1], N)
    basisCoefficients[i-1] = copy(alphas) # copy

    # 3) use exercise strategy to divide paths into exercise and non-exercise domains
    for j in eachindex(exerciseData)
      if exerciseData[j].isValid
        exerciseValue = exerciseData[j].exerciseValue
        continuationValue = exerciseData[j].cumulatedCashFlows
        estimatedContinuationValue = dot(exerciseData[j].values, alphas)

        # for exercise paths, add deflated rebate to
        # deflated cash-flows at previous time frame;
        # for non-exercise paths, add deflated cash-flows to
        # deflated cash-flows at previous time frame
        value = estimatedContinuationValue <= exerciseValue ? exerciseValue : continuationValue

        simulationData[i-1][j].cumulatedCashFlows += value
      end
    end
  end

  # the value of the prodcut can now be estimated by averaging over all paths
  estimatedData = simulationData[1]
  estimate = gen_RiskStatistics(length(estimatedData))

  for j in eachindex(estimatedData)
    add_sample!(estimate, estimatedData[j].cumulatedCashFlows, j)
  end

  return stats_mean(estimate)
end
