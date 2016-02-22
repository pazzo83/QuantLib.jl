type HestonGaussLaguerre <: HestonIntegration end

type Gatheral <: ComplexLogFormula end

type AnalyticHestonEngine{I <: Integer, C <: ComplexLogFormula, HI <: HestonIntegration} <: PricingEngine
  model::HestonModel
  evaluations::I
  cpxLog::C
  integration::HI
end
