struct JamshidianSwaptionEngine{S <: ShortRateModel} <: PricingEngine
  model::S
  # ts::Y

  # function call{S}(::Type{JamshidianSwaptionEngine}, model::S)
  #   new{S, YieldTermStructure}(model)
  # end
  #
  # function call{S, Y}(::Type{JamshidianSwaptionEngine}, model::S, yts::Y)
  #   new{S, Y}(model, yts)
  # end
end

# methods #
function _calculate!(pe::JamshidianSwaptionEngine, swaption::Swaption)
  tsmodel = pe.model
  ref_date = reference_date(tsmodel.ts)
  dc = tsmodel.ts.dc

  amounts = copy(swaption.swap.args.fixedCoupons)
  amounts[end] += swaption.swap.nominal

  maturity = year_fraction(dc, ref_date, swaption.exercise.dates[1])

  fixedPayTimes = zeros(length(swaption.swap.args.fixedPayDates))
  valueTime = year_fraction(dc, ref_date, swaption.swap.args.fixedResetDates[1])
  # for i = 1:length(fixedPayTimes)
  #   fixedPayTimes[i] = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])
  # end
  map!(x -> year_fraction(dc, ref_date, x), fixedPayTimes, swaption.swap.args.fixedPayDates)

  finder = RStarFinder(tsmodel, swaption.swap.nominal, maturity, valueTime, fixedPayTimes, amounts)

  minStrike = -10.0
  maxStrike = 10.0
  slv = BrentSolver(10000, true, true, minStrike, maxStrike)

  rStar = solve(slv, finder, 1e-8, 0.05, minStrike, maxStrike)

  w = isa(swaption.swap.swapT, Payer) ? Put() : Call()

  _size = length(swaption.swap.args.fixedCoupons)

  val = 0.0
  _B = discount_bond(tsmodel, maturity, valueTime, rStar)

  @simd for i = 1:_size
    @inbounds fixedPayTime = year_fraction(dc, ref_date, swaption.swap.args.fixedPayDates[i])

    strike = discount_bond(tsmodel, maturity, fixedPayTime, rStar) / _B

    dboValue = discount_bond_option(tsmodel, w, strike, maturity, valueTime, fixedPayTime)

    @inbounds val += amounts[i] * dboValue
  end

  swaption.results.value = val

  return swaption
end
