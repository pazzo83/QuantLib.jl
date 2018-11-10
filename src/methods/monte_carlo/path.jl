import Base.getindex, Base.setindex!, Base.length, Base.lastindex, Base.copy

struct Path
  tg::TimeGrid
  values::Vector{Float64}
  # TODO ensure values size is the same as tg size
end

Path(tg::TimeGrid) = Path(tg, Vector{Float64}(undef, length(tg.times)))

getindex(p::Path, i::Int) = p.values[i]
setindex!(p::Path, x::Float64, i::Int) = p.values[i] = x
length(p::Path) = length(p.tg.times)
lastindex(p::Path) = lastindex(p.values)

copy(p::Path) = Path(p.tg, copy(p.values))
