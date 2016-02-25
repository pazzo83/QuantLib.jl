type BatesEngine{I <: Integer, C <: ComplexLogFormula, HI <: HestonIntegration} <: AbstractHestonEngine
  model::BatesModel
  evaluations::I
  cpxLog::C
  integration::HI
end

function BatesEngine(batesModel::BatesModel)
  evals = 0
  cpxLog = Gatheral()
  integration = HestonGaussLaguerre(144)

  return BatesEngine(batesModel, evals, cpxLog, integration)
end

function add_on_term(pe::BatesEngine, phi::Float64, t::Float64, j::Int)
  batesModel = pe.model
  nu = get_nu(batesModel)
  delta2 = 0.5 * get_delta(batesModel) * get_delta(batesModel)
  lambda = get_lambda(batesModel)
  i = j == 1 ? 1.0 : 0.0
  g = complex(i, phi)

  return t * lambda * (exp(nu * g + delta2 * g * g) - 1.0 - g * (exp(nu * delta2) - 1.0))
end
