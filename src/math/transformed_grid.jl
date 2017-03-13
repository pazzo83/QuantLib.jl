type TransformedGrid
  grid::Vector{Float64}
  transformedGrid::Vector{Float64}
  dxm::Vector{Float64}
  dxp::Vector{Float64}
  dx::Vector{Float64}
end

function TransformedGrid(grid::Vector{Float64}, f::Function)
  n = length(grid)
  transformedGrid = f(grid)
  dxm = zeros(n)
  dxp = zeros(n)
  dx = zeros(n)

  @inbounds @simd for i = 2:n - 1
    dxm[i] = transformedGrid[i] - transformedGrid[i-1]
    dxp[i] = transformedGrid[i+1] - transformedGrid[i]
    dx[i] = dxm[i] + dxp[i]
  end

  return TransformedGrid(grid, transformedGrid, dxm, dxp, dx)
end

LogGrid(grid::Vector{Float64}) = TransformedGrid(grid, log)
