using StatsFuns

struct FdmMesherComposite{FM1D <: Fdm1DMesher} <: FdmMesher
  layout::FdmLinearOpLayout
  meshers::Vector{FM1D} # this could change
end

# Constructors
function FdmMesherComposite{F1D <: Fdm1DMesher}(mesh::F1D)
  meshers = F1D[mesh]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite{F1D}(layout, meshers)
end

function FdmMesherComposite{F1D <: Fdm1DMesher}(xmesher::F1D, ymesher::F1D)
  meshers = F1D[xmesher, ymesher]
  layout = get_layout_from_meshers(meshers)

  return FdmMesherComposite{F1D}(layout, meshers)
end

function iter_coords!(coord::Vector{Int}, dims::Vector{Int})
  @inbounds for i in eachindex(dims)# = 1:length(dims)
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
  @inbounds @simd for i in eachindex(retVal) # = 1:length(retVal)
    retVal[i] = mesher.meshers[direction].locations[coords[direction]]
    iter_coords!(coords, mesher.layout.dim)
  end

  return retVal
end

get_location(mesher::FdmMesherComposite, coords::Vector{Int}, direction::Int) = mesher.meshers[direction].locations[coords[direction]]

get_dminus(mesher::FdmMesherComposite, coords::Vector{Int}, direction::Int) = mesher.meshers[direction].dminus[coords[direction]]
get_dplus(mesher::FdmMesherComposite, coords::Vector{Int}, direction::Int) = mesher.meshers[direction].dplus[coords[direction]]

## Meshers ##
struct FdmSimpleProcess1dMesher{P <: StochasticProcess1D} <: Fdm1DMesher
  size::Int
  process::P
  maturity::Float64
  tAvgSteps::Int
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
  @inbounds @simd for l = 1:tAvgSteps
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
  @inbounds @simd for i = 1:sz - 1
    dminus[i + 1] = dplus[i] = locations[i + 1] - locations[i]
  end

  dplus[end] = dminus[1] = -1.0

  return FdmSimpleProcess1dMesher(sz, process, maturity, tAvgSteps, _eps, mandatoryPoint, locations, dplus, dminus)
end

struct Concentrating1dMesher <: Fdm1DMesher
  locations::Vector{Float64}
  dplus::Vector{Float64}
  dminus::Vector{Float64}
end

function Concentrating1dMesher(s::Float64, e::Float64, sz::Int64,
  cPoints::Tuple{Float64,Float64}=(typemax(Float64), typemax(Floa64)),
  requireCPoint::Bool=false)
  locations = zeros(sz)
  dplus = zeros(sz)
  dminus = zeros(sz)
  e > s || error("end must be lager than start")

  cPoint = cPoints[1]
  density = cPoints[2]*(e-s)

  ( cPoint == typemax(Float64) || (cPoint >= s && cPoint <= e) ) || error("cPoint must be between start and end")
  ( density == typemax(Float64) || dentity > 0.0 ) || error("density > 0 required")
  ( cPoint == typemax(Float64) || dentity != typemax(Float64) ) || error("density must be given if cPoint is given")
  ( !requireCPoint || cPoint != typemax(Float64) ) || error("cPoint is required in grid but not given")

  dx = 1.0/(sz - 1)
  if cPoint != typemax(Float64)
    u, z = Vector{Float64}(), Vector{Float64}()
    c1 = asinh((s-cPoint)/density)
    c2 = asinh((e-cPoint)/density)
    if requireCPoint
      push!(u, 0.0)
      push!(z, 0.0)
      if !is_close(cPoint, s) && !is_close(cPoint, e)
        z0 = -c1 / (c2 -c1)
        u0 = max( min( (z0 * (sz - 1) + 0.5), (sz) - 2 ), 1) / (size - 1)
        push!(u, u0)
        push!(z, z0)
      end
      push!(u, 1.0)
      push!(z, 1.0)
      transform = LinearInterpolation(u, z)
    end

    li = requireCPoint ? transform.((1:sz)*dx) : (1:sz)*dx
    locations = cPoint + density*sinh(c1*(1.0 - li) + c2*li)
  else
    locations = s + (1:sz)*dx*(e-s)
  end
  locations[1] = locations[end] = 0.0
end
