type TwoFactorShortRateTree{S <: ShortRateDynamics, P1 <: StochasticProcess, P2 <: StochasticProcess} <: ShortRateTree
  tree1::TrinomialTree{P1}
  tree2::TrinomialTree{P2}
  dynamics::S
  treeLattice::TreeLattice2D{TwoFactorShortRateTree{S, P1, P2}}

  function TwoFactorShortRateTree{S, P1, P2}(tree1::TrinomialTree{P1}, tree2::TrinomialTree{P2}, dyn::S)
    twoFactorTree = new{S, P1, P2}(tree1, tree2, dyn)
    twoFactorTree.treeLattice = TreeLattice2D(tree1, tree2, dyn.correlation, twoFactorTree)

    return twoFactorTree
  end
end

function tree(model::TwoFactorModel, grid::TimeGrid)
  #stuff
  dyn = get_dynamics(model)

  tree1 = TrinomialTree(dyn.xProcess, grid)
  tree2 = TrinomialTree(dyn.yProcess, grid)

  return TwoFactorShortRateTree{typeof(dyn), typeof(dyn.xProcess), typeof(dyn.yProcess)}(tree1, tree2, dyn)
end

function rebuild_lattice!(lattice::TwoFactorShortRateTree, tg::TimeGrid)
  rebuild_tree!(lattice.tree1, tg)
  rebuild_tree!(lattice.tree2, tg)

  # update tree lattice
  lattice.treeLattice.tg = tg
  return lattice
end

get_size(tr::TwoFactorShortRateTree, i::Int) = get_size(tr.treeLattice, i)

descendant(tr::TwoFactorShortRateTree, i::Int, idx::Int, branch::Int) = descendant(tr.treeLattice, i, idx, branch)
probability(tr::TwoFactorShortRateTree, i::Int, idx::Int, branch::Int) = probability(tr.treeLattice, i, idx, branch)

function discount(tr::TwoFactorShortRateTree, i::Int, idx::Int)
  modulo = get_size(tr.tree1, i)
  new_idx = idx - 1

  index1 = round(Int, floor(new_idx % modulo)) + 1
  index2 = round(Int, floor(new_idx / modulo)) + 1

  x = get_underlying(tr.tree1, i, index1)
  y = get_underlying(tr.tree2, i, index2)

  r = short_rate(tr.dynamics, tr.treeLattice.tg.times[i], x, y)

  # println("$i $idx $r $x $y $modulo $index1 $index2")

  return exp(-r * tr.treeLattice.tg.dt[i])
end

type G2Dynamics{P <: Parameter} <: ShortRateDynamics
  fitting::P
  xProcess::OrnsteinUhlenbeckProcess
  yProcess::OrnsteinUhlenbeckProcess
  correlation::Float64
end

G2Dynamics{P <: Parameter}(fitting::P, a::Float64, sigma::Float64, b::Float64, eta::Float64, rho::Float64) =
          G2Dynamics{P}(fitting, OrnsteinUhlenbeckProcess(a, sigma), OrnsteinUhlenbeckProcess(b, eta), rho)

short_rate(dyn::G2Dynamics, t::Float64, x::Float64, y::Float64) = dyn.fitting(t) + x + y

type G2{AffineModelType, T <: TermStructure} <: TwoFactorModel{AffineModelType}
  modT::AffineModelType
  a::ConstantParameter
  sigma::ConstantParameter
  b::ConstantParameter
  eta::ConstantParameter
  rho::ConstantParameter
  phi::G2FittingParameter
  ts::T
  privateConstraint::PrivateConstraint{ConstantParameter}
  common::ModelCommon
end

function G2{T <: TermStructure}(ts::T, a::Float64 = 0.1, sigma::Float64 = 0.01, b::Float64 = 0.1, eta::Float64 = 0.01, rho::Float64 = -0.75)
  a_const = ConstantParameter([a], PositiveConstraint())
  sigma_const = ConstantParameter([sigma], PositiveConstraint())
  b_const = ConstantParameter([b], PositiveConstraint())
  eta_const = ConstantParameter([eta], PositiveConstraint())
  rho_const = ConstantParameter([rho], BoundaryConstraint(-1.0, 1.0))

  privateConstraint = PrivateConstraint(ConstantParameter[a_const, sigma_const, b_const, eta_const, rho_const])

  phi = G2FittingParameter(a, sigma, b, eta, rho, ts)

  return G2{AffineModelType, T}(AffineModelType(), a_const, sigma_const, b_const, eta_const, rho_const, phi, ts, privateConstraint, ModelCommon())
end

get_eta(m::G2) = m.eta.data[1]
get_rho(m::G2) = m.rho.data[1]

get_params(m::G2) = Float64[get_a(m), get_sigma(m), get_b(m), get_eta(m), get_rho(m)]

generate_arguments!(m::G2) = m.phi = G2FittingParameter(get_a(m), get_sigma(m), get_b(m), get_eta(m), get_rho(m), m.ts)

get_dynamics(m::G2) = G2Dynamics(m.phi, get_a(m), get_sigma(m), get_b(m), get_eta(m), get_rho(m))

function V(m::G2, t::Float64)
  expat = exp(-get_a(m) * t)
  expbt = exp(-get_b(m) * t)
  cx = get_sigma(m) / get_a(m)
  cy = get_eta(m) / get_b(m)
  valuex = cx * cx * (t + (2.0 * expat - 0.5 * expat * expat - 1.5) / get_a(m))
  valuey = cy * cy * (t + (2.0 * expbt - 0.5 * expbt * expbt - 1.5) / get_b(m))
  value = 2.0 * get_rho(m) * cx * cy * (t + (expat - 1.0) / get_a(m) + (expbt - 1.0) / get_b(m) - (expat * expbt - 1.0) / (get_a(m) + get_b(m)))

  return valuex + valuey + value
end

A(m::G2, t::Float64, T::Float64) = discount(m.ts, T) / discount(m.ts, t) * exp(0.5 * (V(m, T - t) - V(m, T) + V(m, t)))
B(m::G2, x::Float64, t::Float64) = (1.0 - exp(-x * t)) / x

type G2SwaptionPricingFunction <: Function
  a::Float64
  sigma::Float64
  b::Float64
  eta::Float64
  rho::Float64
  w::Float64
  startTime::Float64
  fixedRate::Float64
  sigmax::Float64
  sigmay::Float64
  rhoxy::Float64
  mux::Float64
  muy::Float64
  payTimes::Vector{Float64}
  A::Vector{Float64}
  Ba::Vector{Float64}
  Bb::Vector{Float64}
end

function G2SwaptionPricingFunction(model::G2, w::Float64, startTime::Float64, payTimes::Vector{Float64}, fixedRate::Float64)
  a = get_a(model)
  sigma = get_sigma(model)
  b = get_b(model)
  eta = get_eta(model)
  rho = get_rho(model)

  sigmax = sigma * sqrt(0.5 * (-expm1(-2.0 * a * startTime)) / a)
  sigmay = eta * sqrt(0.5 * (-expm1(-2.0 * b * startTime)) / b)
  rhoxy = rho * eta * sigma * (-expm1(-(a + b) * startTime)) / ((a + b) * sigmax * sigmay)

  temp = sigma * sigma / (a * a)
  mux = -((temp + rho * sigma * eta / (a * b)) * (-expm1(-a * startTime)) - 0.5 * temp * (-expm1(-2.0 * a * startTime)) -
        rho * sigma * eta / (b * (a + b)) * (-expm1(-(b + a) * startTime)))

  temp = eta * eta / (b * b)
  muy = -((temp + rho * sigma * eta / (a * b)) * (-expm1(-b * startTime)) - 0.5 * temp * (-expm1(-2.0 * b * startTime)) -
        rho * sigma * eta / (a * (a + b)) * (-expm1(-(b + a) * startTime)))

  n = length(payTimes)
  A_ = zeros(n)
  Ba = zeros(n)
  Bb = zeros(n)
  for i = 1:n
    A_[i] = A(model, startTime, payTimes[i])
    Ba[i] = B(model, a, payTimes[i] - startTime)
    Bb[i] = B(model, b, payTimes[i] - startTime)
  end

  return G2SwaptionPricingFunction(a, sigma, b, eta, rho, w, startTime, fixedRate, sigmax, sigmay, rhoxy, mux, muy, payTimes, A_, Ba, Bb)
end

function (pricingFunc::G2SwaptionPricingFunction)(x::Float64)
  phi = Normal()
  n = length(pricingFunc.payTimes)

  temp = (x - pricingFunc.mux) / pricingFunc.sigmax
  txy = sqrt(1.0 - pricingFunc.rhoxy * pricingFunc.rhoxy)
  lambda = zeros(n)

  for i = 1:n
    tau = i == 1 ? pricingFunc.payTimes[1] - pricingFunc.startTime : pricingFunc.payTimes[i] - pricingFunc.payTimes[i - 1]
    c = i == n ? (1.0 + pricingFunc.fixedRate * tau) : pricingFunc.fixedRate * tau
    lambda[i] = c * pricingFunc.A[i] * exp(-pricingFunc.Ba[i] * x)
  end

  func = SolvingFunction(lambda, pricingFunc.Bb)
  s1d = BrentSolver(1000)

  yb = solve(s1d, operator(func), 1e-6, 0.00, -100.0, 100.0)

  h1 = (yb - pricingFunc.muy) / (pricingFunc.sigmay * txy) - pricingFunc.rhoxy * (x - pricingFunc.mux) / (pricingFunc.sigmax * txy)
  val = cdf(phi, -pricingFunc.w * h1)

  for i = 1:n
    h2 = h1 + pricingFunc.Bb[i] * pricingFunc.sigmay * sqrt(1.0 - pricingFunc.rhoxy * pricingFunc.rhoxy)
    kappa = -pricingFunc.Bb[i] * (pricingFunc.muy - 0.5 * txy * txy * pricingFunc.sigmay * pricingFunc.sigmay * pricingFunc.Bb[i] +
            pricingFunc.rhoxy * pricingFunc.sigmay * (x - pricingFunc.mux) / pricingFunc.sigmax)

    val -= lambda[i] * exp(kappa) * cdf(phi, -pricingFunc.w * h2)
  end

  blah = exp(-0.5 * temp * temp) * val / (pricingFunc.sigmax * sqrt(2.0 * pi))
  return exp(-0.5 * temp * temp) * val / (pricingFunc.sigmax * sqrt(2.0 * pi))
end

function operator(pricingFunc::G2SwaptionPricingFunction)
  phi = Normal()
  n = length(pricingFunc.payTimes)
  function _inner(x::Float64)
    temp = (x - pricingFunc.mux) / pricingFunc.sigmax
    txy = sqrt(1.0 - pricingFunc.rhoxy * pricingFunc.rhoxy)
    lambda = zeros(n)

    for i = 1:n
      tau = i == 1 ? pricingFunc.payTimes[1] - pricingFunc.startTime : pricingFunc.payTimes[i] - pricingFunc.payTimes[i - 1]
      c = i == n ? (1.0 + pricingFunc.fixedRate * tau) : pricingFunc.fixedRate * tau
      lambda[i] = c * pricingFunc.A[i] * exp(-pricingFunc.Ba[i] * x)
    end

    func = SolvingFunction(lambda, pricingFunc.Bb)
    s1d = BrentSolver(1000)

    yb = solve(s1d, operator(func), 1e-6, 0.00, -100.0, 100.0)

    h1 = (yb - pricingFunc.muy) / (pricingFunc.sigmay * txy) - pricingFunc.rhoxy * (x - pricingFunc.mux) / (pricingFunc.sigmax * txy)
    val = cdf(phi, -pricingFunc.w * h1)

    for i = 1:n
      h2 = h1 + pricingFunc.Bb[i] * pricingFunc.sigmay * sqrt(1.0 - pricingFunc.rhoxy * pricingFunc.rhoxy)
      kappa = -pricingFunc.Bb[i] * (pricingFunc.muy - 0.5 * txy * txy * pricingFunc.sigmay * pricingFunc.sigmay * pricingFunc.Bb[i] +
              pricingFunc.rhoxy * pricingFunc.sigmay * (x - pricingFunc.mux) / pricingFunc.sigmax)

      val -= lambda[i] * exp(kappa) * cdf(phi, -pricingFunc.w * h2)
    end

    blah = exp(-0.5 * temp * temp) * val / (pricingFunc.sigmax * sqrt(2.0 * pi))
    return exp(-0.5 * temp * temp) * val / (pricingFunc.sigmax * sqrt(2.0 * pi))
  end

  return _inner
end

function gen_swaption(model::G2, swaption::Swaption, fixedRate::Float64, range::Float64, intervals::Int)
  settlement = reference_date(model.ts)
  dc = model.ts.dc
  startTime = year_fraction(dc, settlement, swaption.swap.args.floatingResetDates[1])
  w = swaption.swap.payer[2]

  fixedPayTimes = zeros(length(swaption.swap.args.fixedPayDates))
  for i = 1:length(fixedPayTimes)
    fixedPayTimes[i] = year_fraction(dc, settlement, swaption.swap.args.fixedPayDates[i])
  end
  func = G2SwaptionPricingFunction(model, w, startTime, fixedPayTimes, fixedRate)

  upper = func.mux + range * func.sigmax
  lower = func.mux - range * func.sigmax

  integrator = SegmentIntegral(intervals)
  return swaption.swap.nominal * w * discount(model.ts, startTime) * integrator(func, lower, upper)
end

discount_bond(model::G2, t::Float64, T::Float64, x::Float64, y::Float64) = A(model, t, T) * exp(-B(model, get_a(model), (T - t)) * x - B(model, get_b(model), (T - t)) * y)
discount_bond(model::G2, now::Float64, maturity::Float64, factors::Vector{Float64}) = discount_bond(model, now, maturity, factors[1], factors[2])
