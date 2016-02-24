# discretizations
type QuadraticExponentialMartingale <: AbstractHestonDiscretization end

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

get_time(process::HestonProcess, d::Date) = year_fraction(process.riskFreeRate.dc, reference_date(process.riskFreeRate), d)
