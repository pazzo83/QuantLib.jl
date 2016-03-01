type PdeBSM{B <: AbstractBlackScholesProcess}
  process::B
end

diffusion(pde::PdeBSM, t::Float64, x::Float64) = diffusion(pde.process, t, x)
drift(pde::PdeBSM, t::Float64, x::Float64) = drift(pde.process, t, x)
function discount(pde::PdeBSM, t::Float64, ::Float64)
  if abs(t) < 1e-8
    t = 0
  end

  return forward_rate(pde.process.riskFreeRate, t, t, ContinuousCompounding(), NoFrequency()).rate
end

type PdeConstantCoeff
  diffusion::Float64
  drift::Float64
  discount::Float64
end

function PdeConstantCoeff(process::AbstractBlackScholesProcess, t::Float64, x::Float64)
  pde = PdeBSM(process)
  diffusion_ = diffusion(pde, t, x)
  drift_ = drift(pde, t, x)
  discount_ = discount(pde, t, x)

  return PdeConstantCoeff(diffusion_, drift_, discount_)
end

# diffusion(pdecc::PdeConstantCoeff, ::Float64, ::Float64) = pdecc.diffusion
# drift(pdecc::PdeConstantCoeff, ::Float64, ::Float64) = pdecc.drift
# discount

function generate_operator!(pdecc::PdeConstantCoeff, t::Float64, tg::TransformedGrid, L::TridiagonalOperator)
  for i = 2:length(tg.grid) - 1
    sigma = pdecc.diffusion
    nu = pdecc.drift
    r = pdecc.discount

    sigma2 = sigma * sigma

    pd = -(sigma2 / tg.dxm[i] - nu) / tg.dx[i]
    pu = -(sigma2 / tg.dxp[i] + nu) / tg.dx[i]
    pm = sigma2 / (tg.dxm[i] * tg.dxp[i]) + r
    set_mid_row!(L, i, pd, pm, pu)
  end

  return L
end
