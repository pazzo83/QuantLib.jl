type LatticeShortRateModelEngineCommon{T <: ShortRateTree}
  tg::TimeGrid
  lattice::T
end

function update!(eng::LatticeShortRateModelEngine)
  if length(eng.common.tg.times) > 0
    eng.common.lattice = tree(eng.model, eng.common.tg)
  end

  return eng
end
