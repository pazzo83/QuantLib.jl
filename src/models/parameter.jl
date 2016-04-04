# call{P <: Parameter}(param::P, t::Float64) = value(param, t)

type ConstantParameter{C <: Constraint} <: Parameter
  data::Vector{Float64}
  constraint::C
end

get_data(c::ConstantParameter) = c.data
set_params!(c::ConstantParameter, i::Int, val::Float64) = c.data[i] = val
# call(c::ConstantParameter, t::Float64) = value(c, t)

type G2FittingParameter{T <: TermStructure} <: Parameter
  a::Float64
  sigma::Float64
  b::Float64
  eta::Float64
  rho::Float64
  ts::T
end

call(g::G2FittingParameter, t::Float64) = value(g, t)

function value(param::G2FittingParameter, t::Float64)
  forward = forward_rate(param.ts, t, t, ContinuousCompounding(), NoFrequency()).rate
  temp1 = param.sigma * (-expm1(-param.a * t)) / param.a
  temp2 = param.eta * (-expm1(-param.b * t)) / param.b

  val = 0.5 * temp1 * temp1 + 0.5 * temp2 * temp2 + param.rho * temp1 * temp2 + forward

  return val
end

type HullWhiteFittingParameter{T <: TermStructure} <: Parameter
  a::Float64
  sigma::Float64
  ts::T
end

call(h::HullWhiteFittingParameter, t::Float64) = value(h, t)

function value(param::HullWhiteFittingParameter, t::Float64)
  forward = forward_rate(param.ts, t, t, ContinuousCompounding(), NoFrequency()).rate
  temp = param.a < sqrt(eps()) ? param.sigma * t : param.sigma * (-expm1(-param.a * t)) / param.a

  return forward + 0.5 * temp * temp
end

type TermStructureFittingParameter{T <: TermStructure} <: Parameter
  times::Vector{Float64}
  values::Vector{Float64}
  ts::T
end

TermStructureFittingParameter{T <: TermStructure}(ts::T) = TermStructureFittingParameter{T}(zeros(0), zeros(0), ts)

call(tsp::TermStructureFittingParameter, t::Float64) = value(tsp, t)

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

type PiecewiseConstantParameter{C <: Constraint} <: Parameter
  times::Vector{Float64}
  constraint::C

  function PiecewiseConstantParameter(times::Vector{Float64}, constraint::C)
    retTimes = push!(times, 0.0)
    return new{C}(retTimes, constraint)
  end
end

PiecewiseConstantParameter{C <: Constraint}(times::Vector{Float64}, constraint::C) = PiecewiseConstantParameter{C}(times, constraint)

# call(p::PiecewiseConstantParameter, t::Float64) = value(p, t)

set_params!(param::PiecewiseConstantParameter, i::Int, val::Float64) = param.times[i] = val
get_data(param::PiecewiseConstantParameter) = param.times

NullParameter{P <: DataType}(_type::P) = _type([0.0], NoConstraint())

test_params(c::ConstantParameter, params::Vector{Float64}) = QuantLib.Math.test(c.constraint, params)
