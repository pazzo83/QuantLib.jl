call{P <: Parameter}(param::P, t::Float64) = value(param, t)

type ConstantParameter <: Parameter
  data::Vector{Float64}
  constraint::Constraint
end

type G2FittingParameter{T <: TermStructure} <: Parameter
  a::Float64
  sigma::Float64
  b::Float64
  eta::Float64
  rho::Float64
  ts::T
end

function value(param::G2FittingParameter, t::Float64)
  forward = forward_rate(param.ts, t, t, ContinuousCompounding(), NoFrequency()).rate
  temp1 = param.sigma * (1.0 - exp(-param.a * t)) / param.a
  temp2 = param.eta * (1.0 - exp(-param.b * t)) / param.b

  val = 0.5 * temp1 * temp1 + 0.5 * temp2 * temp2 + param.rho * temp1 * temp2 + forward

  return val
end

type HullWhiteFittingParameter{T <: TermStructure} <: Parameter
  a::Float64
  sigma::Float64
  ts::T
end

function value(param::HullWhiteFittingParameter, t::Float64)
  forward = forward_rate(param.ts, t, t, ContinuousCompounding(), NoFrequency()).rate
  temp = param.a < sqrt(eps()) ? param.sigma * t : param.sigma * (1.0 - exp(-param.a * t)) / param.a

  return forward + 0.5 * temp * temp
end

type TermStructureFittingParameter{T <: TermStructure} <: Parameter
  times::Vector{Float64}
  values::Vector{Float64}
  ts::T
end

TermStructureFittingParameter{T <: TermStructure}(ts::T) = TermStructureFittingParameter{T}(zeros(0), zeros(0), ts)

function reset_param_impl!(param::TermStructureFittingParameter)
  param.times = zeros(length(param.times))
  param.values = zeros(length(param.values))

  return param
end

function set_params!(param::TermStructureFittingParameter, tm::Float64, val::Float64)
  push!(param.times, tm)
  push!(param.values, val)

  return param
end

function value(param::TermStructureFittingParameter, t::Float64)
  idx = findfirst(param.times, t)
  return param.values[idx]
end

NullParameter{P <: DataType}(_type::P) = _type([0.0], NoConstraint())

test_params(c::ConstantParameter, params::Vector{Float64}) = JQuantLib.Math.test(c.constraint, params)
