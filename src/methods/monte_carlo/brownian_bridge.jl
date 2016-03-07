type BrownianBridge{I <: Integer}
  size_::I
  t::Vector{Float64}
  sqrtdt::Vector{Float64}
  bridgeIndex::Vector{I}
  leftIndex::Vector{I}
  rightIndex::Vector{I}
  leftWeight::Vector{Float64}
  rightWeight::Vector{Float64}
  stdDev::Vector{Float64}
end

function BrownianBridge(tg::TimeGrid)
  size_ = length(tg.times) - 1
  t = zeros(size_)
  for i in eachindex(t)
    t[i] = tg[i + 1]
  end
  bb = BrownianBridge(size_, t, zeros(size_), zeros(Int, size_), zeros(Int, size_), zeros(Int, size_), zeros(size_), zeros(size_), zeros(size_))
  initialize!(bb)

  return bb
end

function initialize!(bb::BrownianBridge)
  bb.sqrtdt[1] = sqrt(bb.t[1])
  for i = 2:bb.size_
    bb.sqrtdt[i] = sqrt(bb.t[i] - bb.t[i-1])
  end

  # mapVec is used to indicate which points are already constructed.
  # if mapVec[i] is zero, path point i is yet unconstructed
  # mapVec[i] - 1 is the index of the variate that constructs path point # i
  mapVec = zeros(Int, bb.size_)

  # the first point in the construction is the global step
  mapVec[bb.size_] = 1
  # the global step is constructed from the first variate
  bb.bridgeIndex[1] = bb.size_
  # the variance of the global step
  bb.stdDev[1] = sqrt(bb.t[bb.size_])
  # the global step to the last point in time is special
  bb.leftWeight[1] = bb.rightWeight[1] = 0.0

  j = 1
  for i = 1:bb.size_
    # find the next unpopulated entry in the map
    while mapVec[j] != 0
      j += 1
    end
    k = j
    # find the next populated entry from there
    while mapVec[j] == 0
      k += 1
    end
    # l - 1 is now the index of the point to be constructed next
    l = j + ((k - 1 - j) >> 1)
    mapVec[l] = i
    # the i-th Gaussian variate will be used to set point l - 1
    bb.bridgeIndex[i] = l
    bb.leftIndex[i] = j
    bb.rightIndex[i] = k

    if j != 1
      bb.leftWeight[i] = (bb.t[k] - bb.t[l]) / (bb.t[k] - bb.t[j - 1])
      bb.rightWeight[i] = (bb.t[l] - bb.t[j-1]) / (bb.t[k] - bb.t[j - 1])
      bb.stdDev[i] = sqrt(((bb.t[l] - bb.t[j - 1]) * (bb.t[k] - bb.t[l])) / (bb.t[k] - bb.t[j-1]))
    else
      bb.leftWeight[i] = (bb.t[k] - bb.t[l]) / bb.t[k]
      bb.rightWeight[i] = bb.t[l] / bb.t[k]
      bb.stdDev[i] = sqrt(bb.t[l] * (bb.t[k] - bb.t[l]) / bb.t[k])
    end
    j = k + 1
    if j > bb.size_
      # wrap around
      j = 1
    end
  end

  return bb
end
