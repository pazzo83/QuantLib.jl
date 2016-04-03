type BinomialVanillaEngine{P <: AbstractBlackScholesProcess, T <: BinomialTreeType} <: AbstractVanillaEngine
  process::P
  timeSteps::Int
  treeClass::Type{T}
end

function _calculate!(pe::BinomialVanillaEngine, opt::VanillaOption)
  rfdc = pe.process.riskFreeRate.dc
  divdc = pe.process.dividendYield.dc
  voldc = pe.process.blackVolatility.dc
  volcal = pe.process.blackVolatility.calendar

  s0 = state_variable(pe.process).value
  s0 > 0.0 || error("negative or null underlying given")

  maturityDate = opt.exercise.dates[end]
  v = black_vol(pe.process.blackVolatility, maturityDate, s0)
  r = zero_rate(pe.process.riskFreeRate, maturityDate, rfdc, ContinuousCompounding(), NoFrequency()).rate
  q = zero_rate(pe.process.dividendYield, maturityDate, divdc, ContinuousCompounding(), NoFrequency()).rate

  refDate = reference_date(pe.process.riskFreeRate)

  # binomial trees with constant coefficient
  flatRiskFree = FlatForwardTermStructure(refDate, r, rfdc)
  flatDividends = FlatForwardTermStructure(refDate, q, divdc)
  flatVol = BlackConstantVol(refDate, volcal, v, voldc)

  payoff = opt.payoff

  maturity = year_fraction(rfdc, refDate, maturityDate)

  bs = GeneralizedBlackScholesProcess(state_variable(pe.process), flatRiskFree, flatDividends, flatVol)

  grid = TimeGrid(maturity, pe.timeSteps)

  tree = pe.treeClass(bs, maturity, pe.timeSteps, payoff.strike)
  lattice = BlackScholesLattice(tree, r, maturity, pe.timeSteps)

  option = DiscretizedVanillaOption(opt, pe.process, grid)

  initialize!(option, lattice.treeLattice, maturity)
  rollback!(option, grid[3])

  va2 = copy(option.common.values)
  length(va2) == 3 || error("Expect 3 nodes in grid at second step")

  p2d, p2m, p2u = va2

  s2u = get_underlying(lattice, 3, 3) # up price
  s2m = get_underlying(lattice, 3, 2) # middle price
  s2d = get_underlying(lattice, 3, 1) # down (low) price

  # calculate gamma by taking the first derivative of the two deltas
  delta2u = (p2u - p2m) / (s2u - s2m)
  delta2d = (p2m - p2d) / (s2m - s2d)
  gamma = (delta2u - delta2d) / ((s2u - s2d) / 2)

  # rollback to the second-last step, get option vals at this point
  rollback!(option, grid[2])
  va = copy(option.common.values)
  length(va) == 2 || error("Expect 2 nodes in grid at first step")

  p1d, p1u = va

  s1u = get_underlying(lattice, 2, 2) # up (high) price
  s1d = get_underlying(lattice, 2, 1) # down (low) price

  delta = (p1u - p1d) / (s1u - s1d)

  # finally, rollback to t0
  rollback!(option, 0.0)

  p0 = present_value(option)

  opt.results.value = p0
  opt.results.delta = delta
  opt.results.gamma = gamma
  opt.results.theta = black_scholes_theta(pe.process, opt.results.value, opt.results.delta, opt.results.gamma)

  return pe, opt
end
