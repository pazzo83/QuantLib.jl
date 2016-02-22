type HestonModel{CalibratedModelType} <: Model{CalibratedModelType}
  modT::CalibratedModelType
  theta::ConstantParameter
  kappa::ConstantParameter
  sigma::ConstantParameter
  rho::ConstantParameter
  v0::ConstantParameter
  process::HestonProcess
  common::ModelCommon
end

function HestonModel(process::HestonProcess)
  theta = ConstantParameter([process.theta], PositiveConstraint())
  kappa = ConstantParameter([process.kappa], PositiveConstraint())
  sigma = ConstantParameter([process.sigma], PositiveConstraint())
  rho = ConstantParameter([process.rho], BoundaryConstraint(-1.0, 1.0))
  v0 = ConstantParameter([process.v0], PositiveConstraint())

  hm = HestonModel(CalibratedModelType(), theta, kappa, sigma, rho, v0, process, ModelCommon())

  generate_arguments!(hm)

  return hm
end

get_theta(hm::HestonModel) = hm.theta.data[1]
get_kappa(hm::HestonModel) = hm.kappa.data[1]
get_sigma(hm::HestonModel) = hm.sigma.data[1]
get_rho(hm::HestonModel) = hm.rho.data[1]
get_v0(hm::HestonModel) = hm.v0.data[1]

generate_arguments!(hm::HestonModel) =
  hm.process = HestonProcess(hm.process.riskFreeRate, hm.process.dividendYield, hm.process.s0,
              get_theta(hm), get_kappa(hm), get_sigma(hm), get_rho(hm), get_v0(hm))
