import Base.getindex, Base.setindex!, Base.length, Base.endof

type Path
  tg::TimeGrid
  values::Vector{Float64}
  # TODO ensure values size is the same as tg size
end

Path(tg::TimeGrid) = Path(tg, Vector{Float64}(length(tg.times)))

getindex(p::Path, i::Int) = p.values[i]
setindex!(p::Path, x::Float64, i::Int) = p.values[i] = x
length(p::Path) = length(p.tg.times)
endof(p::Path) = endof(p.values)
