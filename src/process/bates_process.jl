type BatesProcess <: AbstractHestonProcess
  s0::Quote
  riskFreeRate::Y1
  dividendYield::Y2
  v0::Float64
  kappa::Float64
  theta::Float64
  sigma::Float64
  rho::Float64
  lambda::Float64
  nu::Float64
  delta::Float64
  m::Float64
  disc::D
  hestonDisc::DH
end
