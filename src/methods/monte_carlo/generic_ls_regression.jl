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
    # println(means)
    # error("DIE")
  end
end
