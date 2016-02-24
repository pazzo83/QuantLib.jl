type BatesModel{CalibratedModelType} <: Model{CalibratedModelType}
  modT::CalibratedModelType
  hestonModel::HestonModel
  nu::ConstantParameter
  delta::ConstantParameter
  lambda::ConstantParameter
  process::BatesProcess
  common::ModelCommon
end

function BatesModel(process::BatesProcess)
  hm = HestonModel(process.hestonProcess)
  nu = ConstantParameter([process.nu], NoConstraint())
  delta = ConstantParameter([process.delta], PositiveConstraint())
  lambda = ConstantParameter([process.lambda], PositiveConstraint())

  bm = BatesModel(CalibratedModelType(), hm, nu, delta, lambda, process, ModelCommon())

  generate_arguments!(bm)

  return bm
end

get_theta(bm::BatesModel) = get_theta(bm.hestonModel)
get_kappa(bm::BatesModel) = get_kappa(bm.hestonModel)
get_sigma(bm::BatesModel) = get_sigma(bm.hestonModel)
get_rho(bm::BatesModel) = get_rho(bm.hestonModel)
get_v0(bm::BatesModel) = get_v0(bm.hestonModel)
get_nu(bm::BatesModel) = bm.nu.data[1]
get_delta(bm::BatesModel) = bm.delta.data[1]
get_lambda(bm::BatesModel) = bm.lambda.data[1]

generate_arguments!(bm::BatesModel) =
  bm.process = BatesProcess(get_risk_free_rate(bm.process), get_dividend_yield(bm.process), get_s0(bm.process),
              get_v0(bm), get_kappa(bm), get_theta(bm), get_sigma(bm), get_rho(bm),
              get_nu(bm), get_delta(bm), get_lambda(bm))
