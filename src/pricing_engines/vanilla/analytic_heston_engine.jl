type HestonGaussLaguerre <: HestonIntegration
  integration::GaussLaguerreIntegration
end

HestonGaussLaguerre(intOrder::Int) = HestonGaussLaguerre(GaussLaguerreIntegration(intOrder))

type Gatheral <: ComplexLogFormula end

type AnalyticHestonEngine{I <: Integer, C <: ComplexLogFormula, HI <: HestonIntegration} <: PricingEngine
  model::HestonModel
  evaluations::I
  cpxLog::C
  integration::HI
end

function AnalyticHestonEngine(hestonModel::HestonModel)
  evals = 0
  cpxLog = Gatheral()
  integration = HestonGaussLaguerre(144)

  return AnalyticHestonEngine(hestonModel, evals, cpxLog, integration)
end
