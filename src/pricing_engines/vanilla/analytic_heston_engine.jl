type HestonGaussLaguerre <: HestonIntegration
  integration::GaussLaguerreIntegration
end

HestonGaussLaguerre(intOrder::Int) = HestonGaussLaguerre(GaussLaguerreIntegration(intOrder))

num_evals(integration::HestonGaussLaguerre) = get_order(integration.integration)

type Gatheral <: ComplexLogFormula end

type AnalyticHestonEngine{C <: ComplexLogFormula} <: AbstractHestonEngine{HestonGaussLaguerre}
  model::HestonModel
  evaluations::Int
  cpxLog::C
  integration::HestonGaussLaguerre
end

function AnalyticHestonEngine(hestonModel::HestonModel)
  evals = 0
  cpxLog = Gatheral()
  integration = HestonGaussLaguerre(144)

  return AnalyticHestonEngine(hestonModel, evals, cpxLog, integration)
end

add_on_term(engine::AnalyticHestonEngine, ::Real, ::Real, ::Int) = complex(0.0)

# Helper #
type FJHelper{C <: ComplexLogFormula, A <: AbstractHestonEngine} <: IntegrationFunction
  j::Int
  kappa::Float64
  theta::Float64
  sigma::Float64
  v0::Float64
  cpxLog::C
  term::Float64
  x::Float64
  sx::Float64
  dd::Float64
  sigma2::Float64
  rsigma::Float64
  t0::Float64
  b::Int
  g_km1::Float64
  engine::A
end

(fjh::FJHelper)(phi::Float64) = fjh(phi, fjh.cpxLog)

function (fjh::FJHelper)(phi::Float64, cpxLog::Gatheral)
  rpsig = fjh.rsigma * phi
  t1 = fjh.t0 + complex(0.0, -rpsig)
  d = sqrt(t1 * t1 - fjh.sigma2 * phi * complex(-phi, fjh.j == 1 ? 1.0 : -1.0))

  ex = exp(-d * fjh.term)
  addOnTerm = add_on_term(fjh.engine, phi, fjh.term, fjh.j)

  if phi != 0.0
    if fjh.sigma > 1e-5
      p = (t1 - d) / (t1 + d)
      g = log((1.0 - p * ex) / (1.0 - p))

      return imag(exp(fjh.v0 * (t1 - d) * (1.0 - ex)/ (fjh.sigma2 * (1.0 - ex * p)) +
                  (fjh.kappa * fjh.theta) / fjh.sigma2 * ((t1 - d) * fjh.term - 2.0 * g) +
                  complex(0.0, phi * (fjh.dd - fjh.sx)) + addOnTerm)) / phi
    else
      #
      td = phi / (2.0 * t1) * complex(-phi, fjh.j == 1 ? 1.0 : -1.0)
      p = td * fjh.sigma2 / (t1 + d)
      g = p * (1.0 - ex)
      return imag(fjh.v0 * td * (1.0 - ex) / (1.0 - p * ex) + (fjh.kappa * fjh.theta) *
                  (td * fjh.term - 2.0 * g / fjh.sigma2) + complex(0.0, phi * (fjh.dd - fjh.sx))
                  + addOnTerm) / phi
    end
  else
    # use l'Hospital's rule to get lim_ {phi -> 0}
    if j == 1
      kmr = fjh.rsigma - fjh.kappa
      if abs(kmr) > 1e-7
        return fjh.dd - fjh.sx + (exp(kmr * fjh.term) * fjh.kappa * fjh.theta - fjh.kappa
                * fjh.theta * (kmr * fjh.term + 1.0)) / (2.0 * kmr * kmr) - fjh.v0 *
                (1.0 - exp(kmr * fjh.term)) / (2.0 * kmr)
      else
        # kappa = rho * sigma
        return fjh.dd - fjh.sx + 0.25 * fjh.kappa * fjh.theta * fjh.term * fjh.term + 0.5 * fjh.v0 * fjh.term
      end
    else
      return fjh.dd - fjh.sx - (exp(-fjh.kappa * fjh.term) * fjh.kappa * fjh.theta + fjh.kappa * fjh.theta *
              (fjh.kappa * fjh.term - 1.0)) / (2.0 * fjh.kappa * fjh.kappa) - fjh.v0 * (1.0 -
              exp(-fjh.kappa * fjh.term)) / (2.0 * fjh.kappa)
    end
  end
end


function FJHelper(kappa::Float64, theta::Float64, sigma::Float64, v0::Float64, s0::Float64, rho::Float64,
          engine::AbstractHestonEngine, cpxLog::ComplexLogFormula, term::Float64, strike::Float64,
          ratio::Float64, j::Int)

  x = log(s0)
  dd = x - log(ratio)
  # sig = get_sigma(engine.model)
  t0 = kappa - (j == 1 ? rho * sigma : 0.0)

  FJHelper(j, kappa, theta, sigma, v0, cpxLog, term, x, log(strike), dd, sigma * sigma, rho * sigma, t0, 0, 0.0, engine)
end

calculate!(integration::HestonGaussLaguerre, ::Float64, f::FJHelper) = integration.integration(f)

function _calculate!(pe::AbstractHestonEngine, opt::EuropeanOption)
  payoff = opt.payoff
  process = pe.model.process

  riskFreeDiscount = discount(get_risk_free_rate(process), opt.exercise.dates[end])
  dividendDiscount = discount(get_dividend_yield(process), opt.exercise.dates[end])

  spotPrice = get_s0(process).value
  strikePrice = payoff.strike

  term = get_time(process, opt.exercise.dates[end])

  val = do_calculation!(pe, opt, riskFreeDiscount, dividendDiscount, spotPrice, strikePrice, term, get_kappa(pe.model),
                  get_theta(pe.model), get_sigma(pe.model), get_v0(pe.model), get_rho(pe.model), payoff,
                  pe.integration, pe.cpxLog, pe.evaluations)

  opt.results.value = val
  return opt
end

get_val(::Call, spotPrice::Float64, dividendDiscount::Float64, riskFreeDiscount::Float64, strikePrice::Float64, p1::Float64, p2::Float64) =
        spotPrice * dividendDiscount * (p1 + 0.5) - strikePrice * riskFreeDiscount * (p2 + 0.5)

get_val(::Put, spotPrice::Float64, dividendDiscount::Float64, riskFreeDiscount::Float64, strikePrice::Float64, p1::Float64, p2::Float64) =
        spotPrice * dividendDiscount * (p1 - 0.5) - strikePrice * riskFreeDiscount * (p2 - 0.5)

function do_calculation!(pe::AbstractHestonEngine, opt::EuropeanOption, riskFreeDiscount::Float64, dividendDiscount::Float64,
                        spotPrice::Float64, strikePrice::Float64, term::Float64, kappa::Float64, theta::Float64,
                        sigma::Float64, v0::Float64, rho::Float64, payoff::PlainVanillaPayoff, integration::HestonGaussLaguerre,
                        cpxLog::ComplexLogFormula, evaluations::Int)
  # calculation stuff
  ratio = riskFreeDiscount / dividendDiscount
  c_inf = min(10.0, max(0.0001, sqrt(1.0 - rho^2) / sigma)) * (v0 + kappa * theta * term)

  evaluations = 0
  p1 = calculate!(integration, c_inf, FJHelper(kappa, theta, sigma, v0, spotPrice, rho, pe, cpxLog, term, strikePrice,
                  ratio, 1)) / pi
  evaluations += num_evals(integration)

  p2 = calculate!(integration, c_inf, FJHelper(kappa, theta, sigma, v0, spotPrice, rho, pe, cpxLog, term, strikePrice,
                  ratio, 2)) / pi
  evaluations += num_evals(integration)

  return get_val(payoff.optionType, spotPrice, dividendDiscount, riskFreeDiscount, strikePrice, p1, p2)
end
