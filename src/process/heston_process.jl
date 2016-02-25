# discretizations
type QuadraticExponentialMartingale <: AbstractHestonDiscretization end
type FullTruncation <: AbstractHestonDiscretization end

type HestonProcess{Y1 <: YieldTermStructure, Y2 <: YieldTermStructure, D <: AbstractDiscretization, DH <: AbstractHestonDiscretization} <: AbstractHestonProcess
  s0::Quote
  riskFreeRate::Y1
  dividendYield::Y2
  v0::Float64
  kappa::Float64
  theta::Float64
  sigma::Float64
  rho::Float64
  disc::D
  hestonDisc::DH
end

function HestonProcess{Y1 <: YieldTermStructure, Y2 <: YieldTermStructure, D <: AbstractHestonDiscretization}(riskFreeRate::Y1,
                      dividendYield::Y2, s0::Quote, v0::Float64, kappa::Float64, theta::Float64, sigma::Float64, rho::Float64,
                      d::D = QuadraticExponentialMartingale())
  disc = EulerDiscretization()
  return HestonProcess{Y1, Y2, EulerDiscretization, D}(s0, riskFreeRate, dividendYield, v0, kappa, theta, sigma, rho, disc, d)
end

get_time(process::AbstractHestonProcess, d::Date) = year_fraction(get_risk_free_rate(process).dc, reference_date(get_risk_free_rate(process)), d)

get_risk_free_rate(hp::HestonProcess) = hp.riskFreeRate
get_dividend_yield(hp::HestonProcess) = hp.dividendYield
get_s0(hp::HestonProcess) = hp.s0
get_v0(hp::HestonProcess) = hp.v0
get_kappa(hp::HestonProcess) = hp.kappa
get_theta(hp::HestonProcess) = hp.theta
get_sigma(hp::HestonProcess) = hp.sigma
get_rho(hp::HestonProcess) = hp.rho
