type BatesProcess{Y1 <: YieldTermStructure, Y2 <: YieldTermStructure, D <: AbstractDiscretization, DH <: AbstractHestonDiscretization} <: AbstractHestonProcess
  hestonProcess::HestonProcess{Y1, Y2, D, DH}
  lambda::Float64
  nu::Float64
  delta::Float64
  m::Float64
end

function BatesProcess{Y1, Y2, DH}(riskFreeRate::Y1, dividendYield::Y2, s0::Quote, v0::Float64,
                      kappa::Float64, theta::Float64, sigma::Float64, rho::Float64, lambda::Float64, nu::Float64,
                      delta::Float64, d::DH = FullTruncation())

  hp = HestonProcess(riskFreeRate, dividendYield, s0, v0, kappa, theta, sigma, rho, d)
  return BatesProcess{Y1, Y2, EulerDiscretization, DH}(hp, lambda, nu, delta, expm1(nu + 0.5 * delta * delta))
end

get_risk_free_rate(bp::BatesProcess) = bp.hestonProcess.riskFreeRate
get_dividend_yield(bp::BatesProcess) = bp.hestonProcess.dividendYield
get_s0(bp::BatesProcess) = bp.hestonProcess.s0
get_v0(bp::BatesProcess) = bp.hestonProcess.v0
get_kappa(bp::BatesProcess) = bp.hestonProcess.kappa
get_theta(bp::BatesProcess) = bp.hestonProcess.theta
get_sigma(bp::BatesProcess) = bp.hestonProcess.sigma
get_rho(bp::BatesProcess) = bp.hestonProcess.rho
get_nu(bp::BatesProcess) = bp.nu
get_lambda(bp::BatesProcess) = bp.lambda
get_delta(bp::BatesProcess) = bp.delta
get_m(bp::BatesProcess) = bp.m
