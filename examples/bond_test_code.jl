# building fixed rate bonds
using JQuantLib

function build_bonds(bond_mats::Vector{Date}, bond_rates::Vector{Float64}, tenor::JQuantLib.Time.TenorPeriod, conv::JQuantLib.Time.BusinessDayConvention,
                    rule::JQuantLib.Time.DateGenerationRule, calendar::JQuantLib.Time.BusinessCalendar, dc::JQuantLib.Time.DayCount, freq::JQuantLib.Time.Frequency, issue_date::Date)
  bonds = Vector{FixedRateBondHelper}(length(bond_mats))
  pricing_engine = DiscountingBondEngine()

  for i =1:length(bond_mats)
    term_date = bond_mats[i]
    # rate = bond_rates[i] / 100
    rate = bond_rates[i]
    sched = JQuantLib.Time.Schedule(issue_date, term_date, tenor, conv, conv, rule, true)
    bond_help = FixedRateBondHelper(Quote(100.0), FixedRateBond(0, 100.0, sched, rate, dc, conv, 100.0, issue_date, calendar, pricing_engine))
    bonds[i] = bond_help
  end

  return bonds
end

function build_depos{P <: Dates.Period, DC <: JQuantLib.Time.DayCount, B <: JQuantLib.Time.BusinessDayConvention, C <: JQuantLib.Time.BusinessCalendar, I <: Integer}(depo_quotes::Vector{Float64}, depo_tenors::Vector{P},
                    dc::DC, conv::B, calendar::C, fixing_days::I)
  depos = Vector{DepositRateHelper}(length(depo_quotes))
  for i = 1:length(depo_quotes)
    depo_quote = Quote(depo_quotes[i])
    depo_tenor = JQuantLib.Time.TenorPeriod(depo_tenors[i])
    depo = DepositRateHelper(depo_quote, depo_tenor, fixing_days, calendar, conv, true, dc)
    depos[i] = depo
  end

  return depos
end

function build_swaps{P <: Dates.Period, DC <: JQuantLib.Time.DayCount, B <: JQuantLib.Time.BusinessDayConvention, C <: JQuantLib.Time.BusinessCalendar, F <: JQuantLib.Time.Frequency, I <: Integer}(swap_quotes::Vector{Float64},
                    swap_tenors::Vector{P}, fixed_dc::DC, fixed_conv::B, calendar::C, fixed_freq::F, float_index::IborIndex, forward_start::I)
  forward_start_period = Dates.Day(forward_start)
  swaps = Vector{SwapRateHelper}(length(swap_quotes))
  pricing_engine = DiscountingSwapEngine()

  for i = 1:length(swap_quotes)
    swaps[i] = SwapRateHelper(swap_quotes[i], swap_tenors[i], calendar, fixed_freq, fixed_conv, fixed_dc, float_index, 0.0, forward_start_period, pricing_engine)
  end

  return swaps
end

function get_npvs(bonds, issue_date, calendar, dc, freq)
  pricing_engine = DiscountingBondEngine()
  rate_quote = Quote(0.05)
  compounding = CompoundedCompounding()
  yts = FlatForwardTermStructure(0, issue_date, calendar, rate_quote, dc, compounding, freq)
  npvs = zeros(length(bonds))
  for i=1:length(bonds)
    npvs[i] = calculate(pricing_engine, yts, bonds[i])
  end

  return npvs
end

function par_rate(yts::YieldTermStructure, dates::Vector{Date}, dc::JQuantLib.Time.DayCount)
  sum = 0.0
  for i = 2:length(dates)
    dt = JQuantLib.Time.year_fraction(dc, dates[i - 1], dates[i])
    sum += discount(yts, dates[i]) * dt
  end

  result = discount(yts, dates[1]) - discount(yts, dates[end])
  return result/sum
end

function setup()
  today = now()
  issue_date = Date(Dates.Year(today), Dates.Month(today), Dates.Day(today))
  bond_mats = [issue_date + Dates.Year(i) for i in range(2, 2, 15)]
  # bond_rates = [5.75, 6.0, 6.25, 6.5, 6.75, 6.80, 7.00, 7.1, 7.15, 7.2]
  bond_rates = [0.0200, 0.0225, 0.0250, 0.0275, 0.0300, 0.0325, 0.0350, 0.0375, 0.0400, 0.0425, 0.0450, 0.0475, 0.0500, 0.0525, 0.0550]
  set_eval_date!(settings, issue_date)

  freq = JQuantLib.Time.Annual()
  tenor = JQuantLib.Time.TenorPeriod(freq)
  conv = JQuantLib.Time.Unadjusted()
  rule = JQuantLib.Time.DateGenerationBackwards()
  calendar = JQuantLib.Time.USGovernmentBondCalendar()
  dc = JQuantLib.Time.SimpleDayCount()
  bonds = build_bonds(bond_mats, bond_rates, tenor, conv, rule, calendar, dc, freq, issue_date)

  return issue_date, bonds, dc, calendar
end

function piecewise_yld_curve()
  issue_date, bonds, dc, calendar = setup()

  interp = JQuantLib.Math.LogInterpolation()
  trait = Discount()
  bootstrap = IterativeBootstrap()

  yts = PiecewiseYieldCurve(issue_date, bonds, dc, interp, trait, 0.00000000001, bootstrap)

  solver = JQuantLib.Math.BrentSolver()
  solver2 = JQuantLib.Math.FiniteDifferenceNewtonSafe()
  calculate!(IterativeBootstrap(), yts, solver2, solver)

  # println(yts.data)

  for bond in bonds
    date_vec = Vector{Date}(length(bond.cashflows.coupons) + 1)
    date_vec[1] = issue_date
    for (i, cf) in enumerate(bond.cashflows.coupons)
      date_vec[i+1] = date(cf)
    end
    par = par_rate(yts, date_vec, dc)
    println(100.0 * par)
  end
end


  # npvs = get_npvs(bonds, issue_date, calendar, dc, freq)
function fitted_bond_curve_exp()
  issue_date, bonds, dc, calendar = setup()

  esf = ExponentialSplinesFitting(true, length(bonds))
  fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, esf, 1e-10, 5000, 1.0)
  initialize!(fitted_curve)
  calculate!(fitted_curve)

  println(fitted_curve.fittingMethod.commons.minimumCostValue)
  println(fitted_curve.fittingMethod.commons.numberOfIterations)
  println(fitted_curve.fittingMethod.commons.guessSolution)

  disc = discount(fitted_curve, 1.0)
  println(disc)
end

function fitted_bond_curve_simp()
  issue_date, bonds, dc, calendar = setup()

  spf = SimplePolynomialFitting(true, 3, length(bonds))
  fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, spf, 1e-10, 5000, 1.0)
  initialize!(fitted_curve)
  calculate!(fitted_curve)

  # println(fitted_curve.fittingMethod.minimumCostValue)
  println(fitted_curve.fittingMethod.commons.numberOfIterations)
  # println(fitted_curve.fittingMethod.guessSolution)

  disc = discount(fitted_curve, 1.0)
  println(disc)
end

function fitted_bond_curve_ns()
  issue_date, bonds, dc, calendar = setup()

  nsf = NelsonSiegelFitting(length(bonds))
  fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, nsf, 1e-10, 5000, 1.0)
  initialize!(fitted_curve)
  calculate!(fitted_curve)

  # println(fitted_curve.fittingMethod.minimumCostValue)
  println(fitted_curve.fittingMethod.commons.numberOfIterations)
  # println(fitted_curve.fittingMethod.guessSolution)

  disc = discount(fitted_curve, 1.0)
  println(disc)
end

function fitted_bond_curve_sven()
  issue_date, bonds, dc, calendar = setup()

  sf = SvenssonFitting(length(bonds))
  fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, sf, 1e-10, 5000, 1.0)
  initialize!(fitted_curve)
  calculate!(fitted_curve)

  # println(fitted_curve.fittingMethod.minimumCostValue)
  println(fitted_curve.fittingMethod.commons.numberOfIterations)
  # println(fitted_curve.fittingMethod.guessSolution)

  disc = discount(fitted_curve, 1.0)
  println(disc)
end

function fitted_bond_curve_cbspline()
  issue_date, bonds, dc, calendar = setup()

  knots = [-30.0, -20.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
  cbsf = CubicBSplinesFitting(true, knots, length(bonds))
  fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, cbsf, 1e-10, 5000, 1.0)
  initialize!(fitted_curve)
  calculate!(fitted_curve)

  println(fitted_curve.fittingMethod.commons.numberOfIterations)
  disc = discount(fitted_curve, 1.0)
  println(disc)
end

function fitted_bond_curve_all()
  tic()
  issue_date, bonds, dc, calendar = setup()

  esf = ExponentialSplinesFitting(true, length(bonds))
  esf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, esf, 1e-10, 5000, 1.0)
  initialize!(esf_fitted_curve)
  calculate!(esf_fitted_curve)

  spf = SimplePolynomialFitting(true, 3, length(bonds))
  spf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, spf, 1e-10, 5000, 1.0)
  initialize!(spf_fitted_curve)
  calculate!(spf_fitted_curve)

  nsf = NelsonSiegelFitting(length(bonds))
  nsf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, nsf, 1e-10, 5000, 1.0)
  initialize!(nsf_fitted_curve)
  calculate!(nsf_fitted_curve)

  knots = [-30.0, -20.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
  cbsf = CubicBSplinesFitting(true, knots, length(bonds))
  cbsf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, cbsf, 1e-10, 5000, 1.0)
  initialize!(cbsf_fitted_curve)
  calculate!(cbsf_fitted_curve)

  sf = SvenssonFitting(length(bonds))
  sf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, sf, 1e-10, 5000, 1.0)
  initialize!(sf_fitted_curve)
  calculate!(sf_fitted_curve)
  toc()

  println("Exponential Splines: $(esf_fitted_curve.fittingMethod.commons.numberOfIterations)")
  println("Simple Polynomial: $(spf_fitted_curve.fittingMethod.commons.numberOfIterations)")
  println("Nelson Siegel: $(nsf_fitted_curve.fittingMethod.commons.numberOfIterations)")
  println("Cubic B-Splines: $(cbsf_fitted_curve.fittingMethod.commons.numberOfIterations)")
  println("Svensson Fitting: $(sf_fitted_curve.fittingMethod.commons.numberOfIterations)")
end

function generate_floatingrate_bond(yts::YieldTermStructure, disc_yts::YieldTermStructure)
  # Floating Rate bond
  settlement_days = 3
  face_amount = 100.0
  fb_issue_date = Date(2005, 10, 21)
  bond_engine = DiscountingBondEngine(disc_yts)
  fb_dc = JQuantLib.Time.Actual360()
  conv = JQuantLib.Time.ModifiedFollowing()
  fb_schedule = JQuantLib.Time.Schedule(Date(2005, 10, 21), Date(2010, 10, 21), JQuantLib.Time.TenorPeriod(JQuantLib.Time.Quarterly()),
                                        JQuantLib.Time.Unadjusted(), JQuantLib.Time.Unadjusted(), JQuantLib.Time.DateGenerationBackwards(), false,
                                        JQuantLib.Time.USNYSECalendar())
  fixing_days = 2
  in_arrears = true
  gearings = ones(length(fb_schedule.dates) - 1)
  spreads = fill(0.001, length(fb_schedule.dates) - 1)
  libor_3m = usd_libor_index(JQuantLib.Time.TenorPeriod(Base.Dates.Month(3)), yts)
  floating_bond = FloatingRateBond(settlement_days, face_amount, fb_schedule, libor_3m, fb_dc, conv, fixing_days, fb_issue_date, bond_engine, in_arrears, 100.0, gearings, spreads)

  return floating_bond
end

function generate_fixedrate_bond(yts::YieldTermStructure)
  settlement_days = 3
  face_amount = 100.0

  fx_schedule = JQuantLib.Time.Schedule(Date(2007, 5, 15), Date(2017, 5, 15), JQuantLib.Time.TenorPeriod(JQuantLib.Time.Semiannual()),
                                        JQuantLib.Time.Unadjusted(), JQuantLib.Time.Unadjusted(), JQuantLib.Time.DateGenerationBackwards(), false,
                                        JQuantLib.Time.USGovernmentBondCalendar())

  pe = DiscountingBondEngine(yts)

  fixedrate_bond = FixedRateBond(settlement_days, face_amount, fx_schedule, 0.045, JQuantLib.Time.ISMAActualActual(), JQuantLib.Time.ModifiedFollowing(),
                                100.0, Date(2007, 5, 15), fx_schedule.cal, pe)

  return fixedrate_bond
end

function generate_discounting_ts(sett::Date)
  settlement_date = sett
  freq = JQuantLib.Time.Semiannual()
  tenor = JQuantLib.Time.TenorPeriod(freq)
  conv = JQuantLib.Time.Unadjusted()
  conv_depo = JQuantLib.Time.ModifiedFollowing()
  rule = JQuantLib.Time.DateGenerationBackwards()
  calendar = JQuantLib.Time.USGovernmentBondCalendar()
  dc_depo = JQuantLib.Time.Actual365()
  dc = JQuantLib.Time.ISDAActualActual()
  dc_bond = JQuantLib.Time.ISMAActualActual()
  fixing_days = 3

  # build depos
  depo_rates = [0.0096, 0.0145, 0.0194]
  depo_tens = [Base.Dates.Month(3), Base.Dates.Month(6), Base.Dates.Month(12)]

  # build bonds
  issue_dates = [Date(2005, 3, 15), Date(2005, 6, 15), Date(2006, 6, 30), Date(2002, 11, 15), Date(1987, 5, 15)]
  mat_dates = [Date(2010, 8, 31), Date(2011, 8, 31), Date(2013, 8, 31), Date(2018, 8, 15), Date(2038, 5, 15)]

  coupon_rates = [0.02375, 0.04625, 0.03125, 0.04000, 0.04500]
  market_quotes = [100.390625, 106.21875, 100.59375, 101.6875, 102.140625]

  insts = Vector{BootstrapHelper}(length(depo_rates) + length(issue_dates))
  for i = 1:length(depo_rates)
    depo_quote = Quote(depo_rates[i])
    depo_tenor = JQuantLib.Time.TenorPeriod(depo_tens[i])
    depo = DepositRateHelper(depo_quote, depo_tenor, fixing_days, calendar, conv_depo, true, dc_depo)
    insts[i] = depo
  end

  pricing_engine = DiscountingBondEngine()

  for i =1:length(coupon_rates)
    term_date = mat_dates[i]
    # rate = bond_rates[i] / 100
    rate = coupon_rates[i]
    issue_date = issue_dates[i]
    market_quote = market_quotes[i]
    sched = JQuantLib.Time.Schedule(issue_date, term_date, tenor, conv, conv, rule, true)
    bond = FixedRateBondHelper(Quote(market_quote), FixedRateBond(3, 100.0, sched, rate, dc_bond, conv, 100.0, issue_date, calendar, pricing_engine))
    insts[i + length(depo_rates)] = bond
  end

  interp = JQuantLib.Math.LogInterpolation()
  trait = Discount()
  bootstrap = IterativeBootstrap()

  yts = PiecewiseYieldCurve(settlement_date, insts, dc, interp, trait, 0.00000000001, bootstrap)

  # solver = JQuantLib.Math.BrentSolver()
  # solver2 = JQuantLib.Math.FiniteDifferenceNewtonSafe()
  calculate!(yts)

  return yts
end

function main()
  issue_date, bonds, dc, calendar = setup()

  println("Today's date: $issue_date")
  println("Calculating fit for 15 bonds....")

  interp = JQuantLib.Math.LogInterpolation()
  trait = Discount()
  bootstrap = IterativeBootstrap()

  yts = PiecewiseYieldCurve(issue_date, bonds, dc, interp, trait, 0.00000000001, bootstrap)

  # solver = JQuantLib.Math.BrentSolver()
  # solver2 = JQuantLib.Math.FiniteDifferenceNewtonSafe()
  calculate!(yts)

  esf = ExponentialSplinesFitting(true, length(bonds))
  esf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, esf, 1e-10, 5000, 1.0)
  initialize!(esf_fitted_curve)
  calculate!(esf_fitted_curve)

  println("(a) exponential splines")
  println("reference date : ", esf_fitted_curve.referenceDate)
  println("number of iterations : ", esf_fitted_curve.fittingMethod.commons.numberOfIterations)

  spf = SimplePolynomialFitting(true, 3, length(bonds))
  spf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, spf, 1e-10, 5000, 1.0)
  initialize!(spf_fitted_curve)
  calculate!(spf_fitted_curve)

  println("(b) simple polynomial")
  println("reference date : ", spf_fitted_curve.referenceDate)
  println("number of iterations : ", spf_fitted_curve.fittingMethod.commons.numberOfIterations)

  nsf = NelsonSiegelFitting(length(bonds))
  nsf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, nsf, 1e-10, 5000, 1.0)
  initialize!(nsf_fitted_curve)
  calculate!(nsf_fitted_curve)

  println("(c) Nelson-Siegel")
  println("reference date : ", nsf_fitted_curve.referenceDate)
  println("number of iterations : ", nsf_fitted_curve.fittingMethod.commons.numberOfIterations)

  knots = [-30.0, -20.0, 0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0]
  cbsf = CubicBSplinesFitting(true, knots, length(bonds))
  cbsf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, cbsf, 1e-10, 5000, 1.0)
  initialize!(cbsf_fitted_curve)
  calculate!(cbsf_fitted_curve)

  println("(d) cubic B-splines")
  println("reference date : ", cbsf_fitted_curve.referenceDate)
  println("number of iterations : ", cbsf_fitted_curve.fittingMethod.commons.numberOfIterations)

  sf = SvenssonFitting(length(bonds))
  sf_fitted_curve = FittedBondDiscountCurve(0, issue_date, calendar, bonds, dc, sf, 1e-10, 5000, 1.0)
  initialize!(sf_fitted_curve)
  calculate!(sf_fitted_curve)

  println("(e) Svensson")
  println("reference date : ", sf_fitted_curve.referenceDate)
  println("number of iterations : ", sf_fitted_curve.fittingMethod.commons.numberOfIterations)

  println("Output par rates for each curve.  In this case, par rates should equal coupons for these par bonds")

  println(" tenor | coupon | bstrap |  (a)  |  (b)  |  (c)  |  (d)  |  (e)  ")

  for bh in bonds
    bond = bh.bond
    date_vec = Vector{Date}(length(bond.cashflows.coupons) + 1)
    date_vec[1] = issue_date
    for (i, cf) in enumerate(bond.cashflows.coupons)
      date_vec[i+1] = date(cf)
    end

    tenor = JQuantLib.Time.year_fraction(dc, issue_date, date(bond.cashflows.coupons[end]))
    println(@sprintf(" %.2f  | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f | %.3f ",
            tenor, 100.0 * bond.cashflows.coupons[end-1].rate.rate, 100.0 * par_rate(yts, date_vec, dc), 100.0 * par_rate(esf_fitted_curve, date_vec, dc), 100.0 * par_rate(spf_fitted_curve, date_vec, dc),
            100.0 * par_rate(nsf_fitted_curve, date_vec, dc), 100.0 * par_rate(cbsf_fitted_curve, date_vec, dc), 100.0 * par_rate(sf_fitted_curve, date_vec, dc)))
  end
end

function main2()
  settlement_date = Date(2008, 9, 18)
  calendar = JQuantLib.Time.USGovernmentBondCalendar()
  set_eval_date!(settings, settlement_date - Base.Dates.Day(3))

  yts = generate_discounting_ts(settlement_date)
  pe = DiscountingBondEngine(yts)

  # build zero coupon bond
  zcb = ZeroCouponBond(3, calendar, 100.0, Date(2013, 8, 15), JQuantLib.Time.Following(), 116.92, Date(2003, 8, 15), pe)

  # println(npv(zcb, pricing_engine, yts))
  # println(clean_price(zcb))
  # println(dirty_price(zcb))
  return npv(zcb, zcb.pricingEngine), clean_price(zcb), dirty_price(zcb)
end

function main3()
  settlement_date = Date(2008, 9, 18)
  set_eval_date!(settings, settlement_date - Dates.Day(3))
  cal = JQuantLib.Time.TargetCalendar()
  dc = JQuantLib.Time.ISDAActualActual()

  # Build deposits
  depo_quotes = [0.043375, 0.031875, 0.0320375, 0.03385, 0.0338125, 0.0335125]
  depo_tenors = [Dates.Week(1), Dates.Month(1), Dates.Month(3), Dates.Month(6), Dates.Month(9), Dates.Year(1)]
  deposit_dc = JQuantLib.Time.Actual360()
  fixing_days = 3
  biz_conv = JQuantLib.Time.ModifiedFollowing()

  depos = build_depos(depo_quotes, depo_tenors, deposit_dc, biz_conv, cal, fixing_days)

  # build swaps
  fixedLegFreq = JQuantLib.Time.Annual()
  fixedLegConv = JQuantLib.Time.Unadjusted()
  fixedLegDC = JQuantLib.Time.EuroThirty360()
  floatingLegIndex = euribor_index(JQuantLib.Time.TenorPeriod(Base.Dates.Month(6)))
  forwardStart = 1

  swap_quotes = [0.0295, 0.0323, 0.0359, 0.0412, 0.0433]
  swap_tenors = [Dates.Year(2), Dates.Year(3), Dates.Year(5), Dates.Year(10), Dates.Year(15)]

  swaps = build_swaps(swap_quotes, swap_tenors, fixedLegDC, fixedLegConv, cal, fixedLegFreq, floatingLegIndex, forwardStart)

  insts = Vector{BootstrapHelper}(length(swap_quotes) + length(depo_quotes))

  insts[1:length(depo_quotes)] = depos
  insts[length(depo_quotes) + 1: end] = swaps

  interp = JQuantLib.Math.LogInterpolation()
  trait = Discount()
  bootstrap = IterativeBootstrap()

  yts = PiecewiseYieldCurve(settlement_date, insts, dc, interp, trait, 1e-15, bootstrap)

  # solver = JQuantLib.Math.BrentSolver()
  # solver2 = JQuantLib.Math.FiniteDifferenceNewtonSafe()
  calculate!(yts)

  disc_yts = generate_discounting_ts(settlement_date)

  # Floating Rate Bond
  fb = generate_floatingrate_bond(yts, disc_yts)

  cap_vol = ConstantOptionVolatility(3, cal, JQuantLib.Time.ModifiedFollowing(), 0.0, JQuantLib.Time.Actual365())
  update_pricer!(fb.cashflows, cap_vol)

  # fb.iborIndex.ts = yts
  # fb.pricingEngine.yts = disc_yts

  # Zero Coupon Bond
  zcb_pe = DiscountingBondEngine(disc_yts)
  zcb_cal = JQuantLib.Time.USGovernmentBondCalendar()

  # build zero coupon bond
  zcb = ZeroCouponBond(3, zcb_cal, 100.0, Date(2013, 8, 15), JQuantLib.Time.Following(), 116.92, Date(2003, 8, 15), zcb_pe)

  # fixed rate bond
  fxb = generate_fixedrate_bond(disc_yts)
  # fxb.pricingEngine.yts = disc_yts

  println("Today's date: ", settings.evaluation_date)
  println("Settlement date: ", settlement_date)
  println("")
  println("              ZC      Fixed    Floating  ")
  println("-----------------------------------------")
  println(@sprintf("  NPV       %.2f   %.2f   %.2f", npv(zcb), npv(fxb), npv(fb)))
  println(@sprintf(" Clean      %.2f   %.2f   %.2f", clean_price(zcb), clean_price(fxb), clean_price(fb)))
  println(@sprintf(" Dirty      %.2f   %.2f   %.2f", dirty_price(zcb), dirty_price(fxb), dirty_price(fb)))
  println(@sprintf("accrued      %.2f     %.2f     %.2f", accrued_amount(zcb, settlement_date), accrued_amount(fxb, settlement_date), accrued_amount(fb, settlement_date)))
  println(@sprintf(" Next C      N/A      %.2f%%    %.2f%%", next_coupon_rate(fxb.cashflows, settlement_date) * 100.0, next_coupon_rate(fb.cashflows, settlement_date) * 100.0))
  println(@sprintf(" Yield      %.2f%%    %.2f%%    %.2f%%",
          JQuantLib.yield(zcb, clean_price(zcb), JQuantLib.Time.Actual360(), CompoundedCompounding(), JQuantLib.Time.Annual(), settlement_date) * 100.0,
          JQuantLib.yield(fxb, clean_price(fxb), JQuantLib.Time.Actual360(), CompoundedCompounding(), JQuantLib.Time.Annual(), settlement_date) * 100.0,
          JQuantLib.yield(fb, clean_price(fb), JQuantLib.Time.Actual360(), CompoundedCompounding(), JQuantLib.Time.Annual(), settlement_date) * 100.0))

  # return npv(fb, fb.pricingEngine), clean_price(fb), dirty_price(fb), JQuantLib.yield(fb, clean_price(fb), JQuantLib.Time.Actual360(), CompoundedCompounding(), JQuantLib.Time.Annual(), settlement_date), next_coupon_rate(fb.cashflows, settlement_date)
end
