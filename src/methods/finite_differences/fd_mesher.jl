using StatsFuns

type FdmMesherComposite{FM1D <: Fdm1DMesher} <: FdmMesher
  layout::FdmLinearOpLayout
  meshers::Vector{FM1D} # this could change
end

# Constructors
function FdmMesherComposite{F1D <: Fdm1DMesher}(mesh::F1D)
  meshers = F1D[mesh]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite(layout, meshers)
end

function FdmMesherComposite{F1D <: Fdm1DMesher}(xmesher::F1D, ymesher::F1D)
  meshers = F1D[xmesher, ymesher]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite(layout, meshers)
end

function iter_coords!(coord::Vector{Int}, dims::Vector{Int})
  for i = 1:length(dims)
    coord[i] += 1
    if coord[i] == dims[i] + 1
      coord[i] = 1
    else
      break
    end
  end

  return coord
end

function get_locations(mesher::FdmMesherComposite, direction::Int)
  coords = ones(Int, length(mesher.layout.dim))
  retVal = zeros(mesher.layout.size)
  for i = 1:length(retVal)
    retVal[i] = mesher.meshers[direction].locations[coords[direction]]
    iter_coords!(coords, mesher.layout.dim)
  end

  return retVal
end

get_location{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].locations[coords[direction]]

get_dminus{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].dminus[coords[direction]]
get_dplus{I <: Integer}(mesher::FdmMesherComposite, coords::Vector{I}, direction::I) = mesher.meshers[direction].dplus[coords[direction]]

## Meshers ##
type FdmSimpleProcess1dMesher{I <: Integer, P <: StochasticProcess1D} <: Fdm1DMesher
  size::I
  process::P
  maturity::Float64
  tAvgSteps::I
  epsilon::Float64
  mandatoryPoint::Float64
  locations::Vector{Float64}
  dplus::Vector{Float64}
  dminus::Vector{Float64}
end

function FdmSimpleProcess1dMesher(sz::Int, process::StochasticProcess1D, maturity::Float64, tAvgSteps::Int, _eps::Float64, mandatoryPoint::Float64 = -1.0)
  locations = zeros(sz)
  dminus = zeros(sz)
  dplus = zeros(sz)
  mp = mandatoryPoint == -1.0 ? process.x0 : mandatoryPoint
  for l = 1:tAvgSteps
    t = (maturity * l) / tAvgSteps

    qMin = min(process.x0, evolve(process, 0.0, process.x0, t, norminvcdf(_eps)))
    qMax = max(max(process.x0, mp), evolve(process, 0.0, process.x0, t, norminvcdf(1.0 - _eps)))

    dp = (1.0 - 2.0 * _eps) / (sz - 1)
    p = _eps
    locations[1] += qMin

    for i = 2:sz - 1
      p += dp
      locations[i] += evolve(process, 0.0, process.x0, t, norminvcdf(p))
    end

    locations[end] += qMax
  end

  locations /= tAvgSteps
  for i = 1:sz - 1
    dminus[i + 1] = dplus[i] = locations[i + 1] - locations[i]
  end

  dplus[end] = dminus[1] = -1.0

  return FdmSimpleProcess1dMesher(sz, process, maturity, tAvgSteps, _eps, mandatoryPoint, locations, dplus, dminus)
end
