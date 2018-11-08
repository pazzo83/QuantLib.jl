struct CachedSwapKey
  index::SwapIndex
  fixing::Date
  tenor::Dates.Period
end

Base.isequal(A::CachedSwapKey, B::CachedSwapKey) = A.index == B.index && A.fixing == B.fixing && A.tenor == B.tenor
# Base.hash(A::CachedSwapKey) = hash(A.index + A.fixing + A.tenor)

mutable struct GSR{TermStructureConsistentModelType, P <: StochasticProcess1D, PARAM1 <: Parameter, PARAM2 <: Parameter, T <: YieldTermStructure} <: Gaussian1DModel{TermStructureConsistentModelType}
  lazyMixin::LazyMixin
  modT::TermStructureConsistentModelType
  stateProcess::P
  evaluationDate::Date
  enforcesTodaysHistoricFixings::Bool
  reversion::PARAM1
  sigma::PARAM2
  volatilities::Vector{Quote}
  reversions::Vector{Quote}
  volstepdates::Vector{Date}
  volsteptimes::Vector{Float64}
  volsteptimesArray::Vector{Float64}
  ts::T
  swapCache::Dict{CachedSwapKey, VanillaSwap}
  common::ModelCommon
end

function GSR(ts::YieldTermStructure, volstepdates::Vector{Date}, volatilities::Vector{Float64}, reversion::Float64, T::Float64 = 60.0)
  volatilityQuotes = Quote[Quote(v) for v in volatilities]
  reversions = [Quote(reversion)]

  volsteptimes, volsteptimesArray = update_times(ts, volstepdates)

  sigma = PiecewiseConstantParameter(volsteptimes, NoConstraint())
  if length(reversions) == 1
    reversion = ConstantParameter([reversions[1].value], NoConstraint())
  else
    reversion = PiecewiseConstantParameter(volsteptimes, NoConstraint())
  end

  @inbounds @simd for i in eachindex(sigma.times)
    set_params!(sigma, i, volatilities[i])
  end

  @inbounds @simd for i in eachindex(get_data(reversion))
    set_params!(reversion, i, reversions[i].value)
  end

  stateProcess = GsrProcess(volsteptimesArray, sigma.times, get_data(reversion), T)

  swapCache = Dict{CachedSwapKey, VanillaSwap}()

  # TOOD registration with quotes in reversions and volatilities

  return GSR(LazyMixin(), TermStructureConsistentModelType(), stateProcess, Date(), true, reversion, sigma, volatilityQuotes, reversions, volstepdates,
        volsteptimes, volsteptimesArray, ts, swapCache, ModelCommon())
end

function update_times(ts::TermStructure, volstepdates::Vector{Date})
  volsteptimes = zeros(length(volstepdates))
  volsteptimesArray = zeros(length(volsteptimes))
  @inbounds @simd for i in eachindex(volstepdates)
    volsteptimes[i] = time_from_reference(ts, volstepdates[i])
    volsteptimesArray[i] = volsteptimes[i]
    if i == 1
      volsteptimes[1] > 0.0 || error("volsteptimes must be positive")
    else
      volsteptimes[i] > volsteptimes[i - 1] || error("volsteptimes must be strictly increasing")
    end
  end

  # TODO flush cache for process

  return volsteptimes, volsteptimesArray
end

function update_times!(model::GSR)
  model.volsteptimes, model.volsteptimesArray = update_times(model.ts, model.volstepdates)
  flush_cache!(model.stateProcess)

  return model
end

function update_state!(model::GSR)
  @simd for i in eachindex(get_data(model.sigma))
    @inbounds set_params!(model.sigma, i, model.volatilities[i].value)
  end
  @simd for i in eachindex(get_data(model.reversion))
    @inbounds set_params!(model.reversion, i, model.reversions[i].value)
  end

  flush_cache!(model.stateProcess)
end

function generate_arguments!(model::GSR)
  flush_cache!(model.stateProcess)
  notify_observers!(model)

  return model
end

function set_params!(model::GSR, params::Vector{Float64})
  if check_params_equal(model, params)
    return model
  end
  n = length(get_data(model.reversion))
  # first set reversion
  @simd for i = 1:n
    @inbounds set_params!(model.reversion, i, params[i])
  end
  # now set sigma
  @simd for i = n + 1:length(get_data(model.sigma))+n
    @inbounds set_params!(model.sigma, i -n, params[i])
  end

  generate_arguments!(model)

  return model
end

# additional methods #
function perform_calculations!(model::GSR)
  if model.evaluationDate == Date() # not calculated
   model.evaluationDate = settings.evaluation_date
  end

  update_times!(model)
  update_state!(model)

  return model
 end

 function get_params(model::GSR)
   # res = Vector{Float64}(length(get_data(model.reversion)) + length(get_data(model.sigma)))
   res = Vector{Float64}()
   append!(res, get_data(model.reversion))
   append!(res, get_data(model.sigma))

   return res
 end

 get_volatilities(model::GSR) = get_data(model.sigma)

 function move_volatility(model::GSR, i::Int)
   res = trues(length(model.reversions) + length(model.volatilities))
   res[length(model.reversions) + i] = false

   return res
 end

 function gaussian_polynomial_integral(model::Gaussian1DModel, a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, y0::Float64, y1::Float64)
   aa = 4.0 * a
   ba = 2.0 * M_SQRT2 * b
   ca = 2.0 * c
   da = M_SQRT2 * d

   x0 = y0 * M_SQRT1_2
   x1 = y1 * M_SQRT1_2

   return (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * erf(x1) - 1.0 / (4.0 * M_SQRTPI) * exp(-x1 * x1) * (2.0 * aa * x1 * x1 * x1 + 3.0 * aa * x1 +
                 2.0 * ba * (x1 * x1 + 1.0) + 2.0 * ca * x1 + 2.0 * da)) - (0.125 * (3.0 * aa + 2.0 * ca + 4.0 * e) * erf(x0) - 1.0 / (4.0 * M_SQRTPI) * exp(-x0 * x0) *
                (2.0 * aa * x0 * x0 * x0 + 3.0 * aa * x0 + 2.0 * ba * (x0 * x0 + 1.0) + 2.0 * ca * x0 + 2.0 * da))
end

function gaussian_shifted_polynomial_integral(model::Gaussian1DModel, a::Float64, b::Float64, c::Float64, d::Float64, e::Float64, h::Float64, x0::Float64, x1::Float64)
   return gaussian_polynomial_integral(model, a, -4.0 * a * h + b, 6.0 * a * h * h - 3.0 * b * h + c, -4 * a * h * h * h + 3.0 * b * h * h - 2.0 * c * h + d,
        a * h * h * h * h - b * h * h * h + c * h * h - d * h + e, x0, x1)
end

function numeraire_impl(model::GSR, t::Float64, y::Float64, yts::YieldTermStructure)
  calculate!(model)

  p = model.stateProcess

  if t == 0.0
    return isa(yts, NullTermStructure) ? discount(model.ts, p.T) : discount(yts, p.T)
  end

  return zerobond(model, p.T, t, y, yts)
end

numeraire(model::Gaussian1DModel, t::Float64, y::Float64, yts::YieldTermStructure) = numeraire_impl(model, t, y, yts)
numeraire(model::Gaussian1DModel, referenceDate::Date, y::Float64, yts::YieldTermStructure) = numeraire(model, time_from_reference(model.ts, referenceDate), y, yts)

function zerobond_impl(model::GSR, T::Float64, t::Float64, y::Float64, yts::YieldTermStructure)
  calculate!(model)

  if t == 0.0
   return isa(yts, NullTermStructure) ? discount(model.ts, T) : discount(yts, T)
  end

  p = model.stateProcess

  # teststDv = std_deviation(p, 0.0, 0.0, t)
  # testEpxt = expectation(p, 0.0, 0.0, t)
  #
  # println(teststDv)
  # println(testEpxt)
  # error("DIE")

  x = y * std_deviation(p, 0.0, 0.0, t) + expectation(p, 0.0, 0.0, t)
  gtT = G!(p, t, T, x)

  d = isa(yts, NullTermStructure) ? discount(model.ts, T) / discount(model.ts, t) : discount(yts, T) / discount(yts, t)

  # println("x: ", x)
  # println("gtT: ", gtT)
  # println("d: ", d)

  return d * exp(-x * gtT - 0.5 * y!(p, t) * gtT * gtT)
end

zerobond(model::Gaussian1DModel, T::Float64, t::Float64, y::Float64, yts::YieldTermStructure) = zerobond_impl(model, T, t, y, yts)
zerobond(model::Gaussian1DModel, maturity::Date, referenceDate::Date, y::Float64, yts::YieldTermStructure) =
        zerobond(model, time_from_reference(yts, maturity), referenceDate != Date() ? time_from_reference(yts, referenceDate) : 0.0, y, yts)

function forward_rate(model::Gaussian1DModel, fixing::Date, referenceDate::Date, y::Float64, iborIdx::IborIndex)
  calculate!(model)

  if fixing <= (model.evaluationDate - Dates.Day(1))
   return fixing(iborIdx, fixing)
  end

  yts = iborIdx.ts

  valueDate = value_date(iborIdx, fixing)
  endDate = advance(iborIdx.tenor.period, iborIdx.fixingCalendar, valueDate, iborIdx.convention)

  dcf = year_fraction(iborIdx.dc, valueDate, endDate)

  return (zerobond(model, valueDate, referenceDate, y, yts) - zerobond(model, endDate, referenceDate, y, yts)) / (dcf * zerobond(model, endDate, referenceDate, y, yts))
end

function y_grid(model::Gaussian1DModel, stdDevs::Float64, gridPoints::Int, T::Float64 = 1.0, t::Float64 = 0.0, y::Float64 = 0.0)
  result = zeros(2 * gridPoints + 1)

  stdDev_0_T = std_deviation(model.stateProcess, 0.0, 0.0, T)
  e_0_T = expectation(model.stateProcess, 0.0, 0.0, T)

  if t < eps()
   stdDev_t_T = stdDev_0_T
   e_t_T = e_0_T
  else
   stdDev_0_t = std_deviation(model.stateProcess, 0.0, 0.0, t)
   stdDev_t_T = std_deviation(model.stateProcess, t, 0.0, T - t)
   e_0_t = expectation(model.stateProcess, 0.0, 0.0, t)
   x_t = y * stdDev_0_t + e_0_t
   e_t_T = expectation(model.stateProcess, t, x_t, T - t)
  end

  h = stdDevs / gridPoints

  @simd for j = -gridPoints:gridPoints
   @inbounds result[j + gridPoints + 1] = (e_t_T + stdDev_t_T * j * h - e_0_T) / stdDev_0_T
  end

  return result
end

function calibrate_volatilities_iterative!(model::GSR, helpers::Vector{H}, method::OptimizationMethod, endCriteria::EndCriteria,
                                          constraint::Constraint = PositiveConstraint(), weights::Vector{Float64} = Vector{Float64}()) where {H <: CalibrationHelper}
  @inbounds @simd for i in eachindex(helpers)
    h = H[helpers[i]]
    calibrate!(model, h, method, endCriteria, constraint, weights, move_volatility(model, i))
  end

  return model
end
