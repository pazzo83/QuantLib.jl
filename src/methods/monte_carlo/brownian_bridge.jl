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

function BrownianBridge(steps::Int)
  t = Vector{Float64}(steps)
  for i in eachindex(t)
    t[i] = float(i)
  end

  bb = BrownianBridge(steps, t, Vector{Float64}(steps), Vector{Int}(steps), Vector{Int}(steps), Vector{Int}(steps), Vector{Float64}(steps), Vector{Float64}(steps), Vector{Float64}(steps))
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
  bb.bridgeIndex[1] = bb.size_ - 1
  # the variance of the global step
  bb.stdDev[1] = sqrt(bb.t[bb.size_])
  # the global step to the last point in time is special
  bb.leftWeight[1] = bb.rightWeight[1] = 0.0

  j = 0
  for i = 1:bb.size_ - 1
    # find the next unpopulated entry in the map
    while mapVec[j+1] != 0
      j += 1
    end
    k = j
    # find the next populated entry from there
    while mapVec[k+1] == 0
      k += 1
    end
    # l - 1 is now the index of the point to be constructed next
    l = j + ((k - 1 - j) >> 1)
    mapVec[l+1] = i
    # the i-th Gaussian variate will be used to set point l - 1
    bb.bridgeIndex[i+1] = l
    bb.leftIndex[i+1] = j
    bb.rightIndex[i+1] = k

    if j != 0
      bb.leftWeight[i+1] = (bb.t[k+1] - bb.t[l+1]) / (bb.t[k+1] - bb.t[j])
      bb.rightWeight[i+1] = (bb.t[l+1] - bb.t[j]) / (bb.t[k+1] - bb.t[j])
      bb.stdDev[i+1] = sqrt(((bb.t[l+1] - bb.t[j]) * (bb.t[k+1] - bb.t[l+1])) / (bb.t[k+1] - bb.t[j]))
    else
      bb.leftWeight[i+1] = (bb.t[k+1] - bb.t[l+1]) / bb.t[k+1]
      bb.rightWeight[i+1] = bb.t[l+1] / bb.t[k+1]
      bb.stdDev[i+1] = sqrt(bb.t[l+1] * (bb.t[k+1] - bb.t[l+1]) / bb.t[k+1])
    end
    j = k + 1
    if j >= bb.size_
      # wrap around
      j = 0
    end
  end

  return bb
end
