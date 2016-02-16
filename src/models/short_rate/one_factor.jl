## Consts (might move these to MATH)
const M_SQRT2 = 1.41421356237309504880168872420969808   # sqrt(2)
const M_SQRT1_2 = 0.707106781186547524400844362104849039  # 1/sqrt(2)
const M_SQRTPI = 1.77245385090551602792981

## ONE FACTOR MODELS ##
type OneFactorShortRateTree{S <: ShortRateDynamics} <: ShortRateTree
  tree::TrinomialTree
  dynamics::S
  tg::TimeGrid
  treeLattice::TreeLattice1D

  function OneFactorShortRateTree{S}(tree::TrinomialTree, dynamics::S, tg::TimeGrid)
    oneFactorTree = new(tree, dynamics, tg)
    oneFactorTree.treeLattice = TreeLattice1D(tg, get_size(tree, 2), oneFactorTree)

    return oneFactorTree
  end
end

get_size{I <: Integer}(tr::OneFactorShortRateTree, i::I) = get_size(tr.tree, i)

function discount(tr::OneFactorShortRateTree, i::Int, idx::Int)
  x = get_underlying(tr.tree, i, idx)
  r = short_rate(tr.dynamics, tr.tg.times[i], x)
  return exp(-r * tr.tg.dt[i])
end

descendant(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = descendant(tr.tree, i, idx, branch)
probability(tr::OneFactorShortRateTree, i::Int, idx::Int, branch::Int) = probability(tr.tree, i, idx, branch)

get_params(m::OneFactorModel) = Float64[get_a(m), get_sigma(m)]

type BlackKarasinski{TermStructureConsistentModelType, T <: TermStructure} <: OneFactorModel{TermStructureConsistentModelType}
  modT::TermStructureConsistentModelType
  a::ConstantParameter
  sigma::ConstantParameter
  ts::T
  privateConstraint::PrivateConstraint
  common::ModelCommon
end

function BlackKarasinski(ts::TermStructure, a::Float64 = 0.1, sigma = 0.1)
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  return BlackKarasinski(TermStructureConsistentModelType(), a_const, sigma_const, ts, privateConstraint, ModelCommon()) # ShortRateModelCommon())
end

generate_arguments!(m::BlackKarasinski) = m # do nothing


type HullWhite{AffineModelType, T <: TermStructure} <: OneFactorModel{AffineModelType}
  modT::AffineModelType
  r0::Float64
  a::ConstantParameter
  sigma::ConstantParameter
  phi::HullWhiteFittingParameter
  ts::T
  privateConstraint::PrivateConstraint
  common::ModelCommon
end

function HullWhite{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01)
  _rate = forward_rate(ts, 0.0, 0.0, ContinuousCompounding(), NoFrequency())
  r0 = _rate.rate
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const])

  phi  = HullWhiteFittingParameter(a, sigma, ts)

  return HullWhite(AffineModelType(), r0, a_const, sigma_const, phi, ts, privateConstraint, ModelCommon()) #, ShortRateModelCommon())
end

immutable CachedSwapKey
  index::SwapIndex
  fixing::Date
  tenor::Dates.Period
end

Base.isequal(A::CachedSwapKey, B::CachedSwapKey) = A.index == B.index && A.fixing == B.fixing && A.tenor == B.tenor
# Base.hash(A::CachedSwapKey) = hash(A.index + A.fixing + A.tenor)

type GSR{TermStructureConsistentModelType, P <: StochasticProcess1D, PARAM1 <: Parameter, PARAM2 <: Parameter, T <: YieldTermStructure} <: Gaussian1DModel{TermStructureConsistentModelType}
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

  for i in eachindex(sigma.times)
    set_params!(sigma, i, volatilities[i])
  end

  for i in eachindex(get_data(reversion))
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
  for i in eachindex(volstepdates)
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
  for i in eachindex(get_data(model.sigma))
    set_params!(model.sigma, i, model.volatilities[i].value)
  end
  for i in eachindex(get_data(model.reversion))
    set_params!(model.reversion, i, model.reversions[i].value)
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
  for i = 1:n
    set_params!(model.reversion, i, params[i])
  end
  # now set sigma
  for i = n + 1:length(get_data(model.sigma))+n
    set_params!(model.sigma, i -n, params[i])
  end

  generate_arguments!(model)

  return model
end

## Dynamics ##

type HullWhiteDynamics{P <: Parameter} <: ShortRateDynamics
  process::OrnsteinUhlenbeckProcess
  fitting::P
  a::Float64
  sigma::Float64
end

HullWhiteDynamics{P <: Parameter}(fitting::P, a::Float64, sigma::Float64) = HullWhiteDynamics(OrnsteinUhlenbeckProcess(a, sigma), fitting, a, sigma)

short_rate(dynamic::HullWhiteDynamics, t::Float64, x::Float64) = x + dynamic.fitting(t)

generate_arguments!(m::HullWhite) = m.phi = HullWhiteFittingParameter(get_a(m), get_sigma(m), m.ts)

get_dynamics(m::HullWhite) = HullWhiteDynamics(m.phi, get_a(m), get_sigma(m))

type BlackKarasinskiDynamics{P <: Parameter} <: ShortRateDynamics
  process::OrnsteinUhlenbeckProcess
  fitting::P
  a::Float64
  sigma::Float64
end

BlackKarasinskiDynamics{P <: Parameter}(fitting::P, a::Float64, sigma::Float64) = BlackKarasinskiDynamics(OrnsteinUhlenbeckProcess(a, sigma), fitting, a, sigma)

short_rate(dynamic::BlackKarasinskiDynamics, t::Float64, x::Float64) = exp(x + dynamic.fitting(t))

immutable BlackKarasinskiHelper{I <: Integer}
  sz::I
  xMin::Float64
  dx::Float64
  dt::Float64
  discountBondPrice::Float64
  statePrices::Vector{Float64}
end

function BlackKarasinskiHelper{I <: Integer}(i::I, xMin::Float64, dx::Float64, discountBond::Float64, _tree::OneFactorShortRateTree)
  sz = get_size(_tree, i)
  statePrices = get_state_prices!(_tree, i)
  dt = _tree.tg.dt[i]

  return BlackKarasinskiHelper(sz, xMin, dx, dt, discountBond, statePrices)
end

function operator(bkh::BlackKarasinskiHelper)
  function _inner(theta::Float64)
    val = bkh.discountBondPrice
    x = bkh.xMin

    for j = 1:bkh.sz
      disc = exp(-exp(theta + x) * bkh.dt)
      val -= bkh.statePrices[j] * disc
      x += bkh.dx
    end

    return val
  end

  return _inner
end

## Tree methods ##

function tree(model::HullWhite, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)
  numericDynamics = HullWhiteDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{HullWhiteDynamics}(trinomial, numericDynamics, grid)

  reset_param_impl!(phi)

  @simd for i = 1:length(grid.times) - 1
    @inbounds discountBond = discount(model.ts, grid.times[i + 1])
    statePrices = get_state_prices!(numericTree, i)
    sz = get_size(numericTree, i)
    @inbounds dt = grid.dt[i]
    @inbounds dx = trinomial.dx[i]
    x = get_underlying(trinomial, i, 1)
    val = 0.0
    for j = 1:sz
      @inbounds val += statePrices[j] * exp(-x * dt)
      x += dx
    end
    val = log(val / discountBond) / dt
    @inbounds set_params!(phi, grid.times[i], val)
  end

  return numericTree
end

function tree(model::BlackKarasinski, grid::TimeGrid)
  phi = TermStructureFittingParameter(model.ts)

  numericDynamics = BlackKarasinskiDynamics(phi, get_a(model), get_sigma(model))
  trinomial = TrinomialTree(numericDynamics.process, grid)
  numericTree = OneFactorShortRateTree{BlackKarasinskiDynamics}(trinomial, numericDynamics, grid)

  reset_param_impl!(phi)

  val = 1.0
  vMin = -50.0
  vMax = 50.0

  @simd for i = 1:length(grid.times) - 1
    @inbounds discountBond = discount(model.ts, grid.times[i + 1])
    xMin = get_underlying(trinomial, i, 1)
    @inbounds dx = trinomial.dx[i]

    solverHelper = BlackKarasinskiHelper(i, xMin, dx, discountBond, numericTree)
    slvr = BrentSolver(1000) # max evals = 1000
    val = solve(slvr, operator(solverHelper), 1e-7, val, vMin, vMax)

    @inbounds set_params!(phi, grid.times[i], val)
  end

  return numericTree
end


type RStarFinder{M <: ShortRateModel}
  model::M
  strike::Float64
  maturity::Float64
  valueTime::Float64
  fixedPayTimes::Vector{Float64}
  amounts::Vector{Float64}
end

function operator(rsf::RStarFinder)
  function _inner(x::Float64)
    _value = rsf.strike
    _B = discount_bond(rsf.model, rsf.maturity, rsf.valueTime, x)
    sz = length(rsf.fixedPayTimes)
    for i = 1:sz
      dbVal = discount_bond(rsf.model, rsf.maturity, rsf.fixedPayTimes[i], x) / _B
      _value -= rsf.amounts[i] * dbVal
    end

    return _value
  end

  return _inner
end

function B(model::HullWhite, t::Float64, T::Float64)
  _a = get_a(model)
  if _a < sqrt(eps())
    return T - t
  else
    return (1.0 - exp(-_a * (T - t))) / _a
  end
end

function A(model::HullWhite, t::Float64, T::Float64)
  discount1 = discount(model.ts, t)
  discount2 = discount(model.ts, T)

  forward = forward_rate(model.ts, t, t, ContinuousCompounding(), NoFrequency())
  temp = get_sigma(model) * B(model, t, T)
  val = B(model, t, T) * forward.rate - 0.25 * temp * temp * B(model, 0.0, 2.0* t)

  return exp(val) * discount2 / discount1
end

# type HullWhite{T <: TermStructure} <: ShortRateModel
#   r0::Float64
#   a::ConstantParameter
#   sigma::ConstantParameter
#   phi::HullWhiteFittingParameter
# end

discount_bond(model::OneFactorModel, tNow::Float64, maturity::Float64, factors::Vector{Float64}) = discount_bond(model, tNow, maturity, factors[1])
discount_bond{M <: OneFactorModel}(model::M, tNow::Float64, maturity::Float64, _rate::Float64) = A(model, tNow, maturity) * exp(-B(model, tNow, maturity) * _rate)

function discount_bond_option{O <: OptionType}(model::HullWhite, optionType::O, strike::Float64, maturity::Float64, bondStart::Float64, bondMaturity::Float64)
  _a = get_a(model)
  if _a < sqrt(eps())
    v = get_sigma(model) * B(model, bondStart, bondMaturity) * sqrt(maturity)
  else
    v = get_sigma(model) / (_a * sqrt(2.0 * _a)) * sqrt(exp(-2.0 * _a * (bondStart - maturity)) - exp(-2.0 * _a * bondStart) - 2.0 * (exp(-_a * (bondStart + bondMaturity - 2.0 * maturity))
        - exp(-_a * (bondStart + bondMaturity))) + exp(-2.0 * _a * (bondMaturity - maturity)) - exp(-2.0 * _a * bondMaturity))
  end

   f = discount(model.ts, bondMaturity)
   k = discount(model.ts, bondStart) * strike

   return black_formula(optionType, k, f, v)
end

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
   e_t_T = expecation(model.stateProcess, t, x_t, T - t)
  end

  h = stdDevs / gridPoints

  for j = -gridPoints:gridPoints
   result[j + gridPoints + 1] = (e_t_T + stdDev_t_T * j * h - e_0_T) / stdDev_0_T
  end

  return result
end

function calibrate_volatilities_iterative!{H <: CalibrationHelper}(model::GSR, helpers::Vector{H}, method::OptimizationMethod, endCriteria::EndCriteria,
                                          constraint::Constraint = PositiveConstraint(), weights::Vector{Float64} = Vector{Float64}())
  for i in eachindex(helpers)
    h = H[helpers[i]]
    calibrate!(model, h, method, endCriteria, constraint, weights, move_volatility(model, i))
  end

  return model
end
