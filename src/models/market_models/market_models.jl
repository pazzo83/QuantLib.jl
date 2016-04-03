## General Market Model methods ##
function get_covariance(mm::AbstractMarketModel, i::Int)
  if isempty(mm.covariance)
    resize!(mm.covariance, mm.numberOfSteps)
    for j = 1:mm.numberOfSteps
      mm.covariance[j] = mm.pseudoRoots[j] * mm.pseudoRoots[j]'
    end
    i <= length(mm.covariance) || error("i must be less than or equal to covariance size")
  end

  return mm.covariance[i]
end

function total_covariance(mm::AbstractMarketModel, endIndex::Int)
  if isempty(mm.totalCovariance)
    resize!(mm.totalCovariance, mm.numberOfSteps)
    mm.totalCovariance[1] = get_covariance(mm, 1) # need to trigger calculation
    for j = 2:mm.numberOfSteps
      mm.totalCovariance[j] = mm.totalCovariance[j-1] + mm.covariance[j]
    end
    endIndex <= length(mm.covariance) || error("endIndex must be less than or equal to total covariance size")
  end

  return mm.totalCovariance[endIndex]
end
