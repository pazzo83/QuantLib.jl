Getting Started
===============

Let's look at a few examples!

First, start off by importing QuantLib.jl

.. code-block:: julia

    using QuantLib



Price a fixed rate bond
-----------------------

First, set up the global environment.

.. code-block:: julia

    settlement_date = Date(2008, 9, 18) # construct settlement date
    # settings is a global singleton that contains global settings
    set_eval_date!(settings, settlement_date - Dates.Day(3))

Then, construct settings we will need to construct the yield curve

.. code-block:: julia

    freq = QuantLib.Time.Semiannual()
    tenor = QuantLib.Time.TenorPeriod(freq)
    conv = QuantLib.Time.Unadjusted()
    conv_depo = QuantLib.Time.ModifiedFollowing()
    rule = QuantLib.Time.DateGenerationBackwards()
    calendar = QuantLib.Time.USGovernmentBondCalendar()
    dc_depo = QuantLib.Time.Actual365()
    dc = QuantLib.Time.ISDAActualActual()
    dc_bond = QuantLib.Time.ISMAActualActual()
    fixing_days = 3

Build the Deposit and Bond helpers we will use to bootstrap the curve.

First we will set up vectors of the deposit rates and tenors, and bond issue dates and maturity dates (these could come from some other market data source)

.. code-block:: julia

    # build depos
    depo_rates = [0.0096, 0.0145, 0.0194]
    depo_tens = [Base.Dates.Month(3), Base.Dates.Month(6), Base.Dates.Month(12)]

    # build bonds
    issue_dates = [Date(2005, 3, 15), Date(2005, 6, 15), Date(2006, 6, 30), Date(2002, 11, 15), Date(1987, 5, 15)]
    mat_dates = [Date(2010, 8, 31), Date(2011, 8, 31), Date(2013, 8, 31), Date(2018, 8, 15), Date(2038, 5, 15)]

Now, we'll use some static coupon rates and market quotes for the bond helpers

.. code-block:: julia

    coupon_rates = [0.02375, 0.04625, 0.03125, 0.04000, 0.04500]
    market_quotes = [100.390625, 106.21875, 100.59375, 101.6875, 102.140625]

With this data, we can now construct our helpers

.. code-block:: julia

    # construct the deposit and fixed rate bond helpers
    insts = Vector{BootstrapHelper}(length(depo_rates) + length(issue_dates))
    for i = eachindex(depo_rates)
      depo_quote = Quote(depo_rates[i])
      depo_tenor = QuantLib.Time.TenorPeriod(depo_tens[i])
      depo = DepositRateHelper(depo_quote, depo_tenor, fixing_days, calendar, conv_depo, true, dc_depo)
      insts[i] = depo
    end

    for i = eachindex(coupon_rates)
      term_date = mat_dates[i]
      rate = coupon_rates[i]
      issue_date = issue_dates[i]
      market_quote = market_quotes[i]
      sched = QuantLib.Time.Schedule(issue_date, term_date, tenor, conv, conv, rule, true)
      bond = FixedRateBondHelper(Quote(market_quote), FixedRateBond(3, 100.0, sched, rate, dc_bond, conv,
                                100.0, issue_date, calendar, DiscountingBondEngine()))
      insts[i + length(depo_rates)] = bond
    end

With our helpers created, we can start to construct the yield curve which we will bootstrap

.. code-block:: julia

    interp = QuantLib.Math.LogInterpolation()
    trait = Discount()
    bootstrap = IterativeBootstrap()
    yts = PiecewiseYieldCurve(settlement_date, insts, dc, interp, trait, 0.00000000001, bootstrap)

Now, we can trigger the bootstrapping calculation (this can be triggered by a number of events, but for now we will just directly trigger calculation)

.. code-block:: julia

    calculate!(yts)

Let's now create our fixed rate bond, by generating a coupon schedule and giving it a pricing engine

.. code-block:: julia

    settlement_days = 3
    face_amount = 100.0

    fixed_schedule = QuantLib.Time.Schedule(Date(2007, 5, 15), Date(2017, 5, 15),
                QuantLib.Time.TenorPeriod(QuantLib.Time.Semiannual()), QuantLib.Time.Unadjusted(),
                QuantLib.Time.Unadjusted(), QuantLib.Time.DateGenerationBackwards(), false,
                QuantLib.Time.USGovernmentBondCalendar())

    pe = DiscountingBondEngine(yts)

    fixedrate_bond = FixedRateBond(settlement_days, face_amount, fixed_schedule, 0.045,
                  QuantLib.Time.ISMAActualActual(), QuantLib.Time.ModifiedFollowing(), 100.0,
                  Date(2007, 5, 15), fixed_schedule.cal, pe)

Finally, we can request for the bond's NPV!

.. code-block:: julia

    npv(fixedrate_bond) # 107.66828913260542



Calculate the Survival Probability of a Credit Default Swap
-----------------------------------------------------------

First, let's set up the environment

.. code-block:: julia

    cal = QuantLib.Time.TargetCalendar()
    todays_date = Date(2007, 5, 15)
    settlementDate = todays_date
    set_eval_date!(settings, todays_date)

Now, let's generate a flat-forward term structure for use with our CDS Helpers (which are used to generate the Hazard Rate Curve)

.. code-block:: julia

    flatRate = Quote(0.01)

    tsCurve = FlatForwardTermStructure(settlementDate, cal, flatRate, QuantLib.Time.Actual365())

To bootstrap the hazard rate curve that we will use for survival probability (and inversely, default probability), we need to build CDS helpers.  To begin, we'll set a recovery rate, and quote spreads, tenors, and maturity dates for 4 CDS helpers

.. code-block:: julia

    recoveryRate = 0.5
    quoteSpreads = [0.0150, 0.0150, 0.0150, 0.0150]
    tenors = [Dates.Month(3), Dates.Month(6), Dates.Year(1), Dates.Year(2)]

    maturities = [QuantLib.Time.adjust(cal, QuantLib.Following(), todays_date + ten) for ten in tenors]

Let's build our CDS helpers

.. code-block:: julia

    insts = SpreadCDSHelper[SpreadCDSHelper(Quote(quoteSpreads[i]), tenors[i], 0, cal, QuantLib.Time.Quarterly(),
            QuantLib.Time.Following(), QuantLib.Time.DateGenerationTwentieth(), QuantLib.Time.Actual365(),
            recoveryRate, tsCurve) for i in eachindex(tenors)]

With our helpers constructed, now we can build the hazard rate curve.

.. code-block:: julia

    hazardRateStructure = PiecewiseDefaultCurve(todays_date, insts, QuantLib.Time.Actual365(),
                          QuantLib.Math.BackwardFlatInterpolation(), HazardRate(), 1.0e-12)

By requested for the curve nodes, we will trigger the bootstrap calculation

.. code-block:: julia

    hr_curve_data = nodes(hazardRateStructure)

Now we can output the 1Y and 2Y survival probabilities

.. code-block:: julia

    println(@sprintf("1Y Survival Probability: %.6f %%", survival_probability(hazardRateStructure,
            todays_date + Dates.Year(1)) * 100.0))
    println(@sprintf("2Y Survival Probability: %.6f %%", survival_probability(hazardRateStructure,
            todays_date + Dates.Year(2)) * 100.0))



Price a Swaption Using a G2 Calibrated Model
-----------------------------------------

Set up our environment

.. code-block:: julia

    cal = QuantLib.Time.TargetCalendar()
    settlementDate = Date(2002, 2, 19)
    todays_date = Date(2002, 2, 15)
    set_eval_date!(settings, todays_date)

Gather appropriate market data

.. code-block:: julia

    swaptionMats = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]
    swaptionVols = [0.1490, 0.1340, 0.1228, 0.1189, 0.1148, 0.1290, 0.1201, 0.1146, 0.1108,
                    0.1040, 0.1149, 0.1112, 0.1070, 0.1010, 0.0957, 0.1047, 0.1021, 0.0980, 0.0951,
                    0.1270, 0.1000, 0.0950, 0.0900, 0.1230, 0.1160]
    swaptionLengths = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]

Generate a flat-forward term structure implying a 1x5 swap at 5%

.. code-block:: julia

    flat_rate = Quote(0.04875825)
    rhTermStructure = FlatForwardTermStructure(settlementDate, cal, flat_rate, QuantLib.Time.Actual365())

Build an ATM swap

.. code-block:: julia

    fixedLegFrequency = QuantLib.Time.Annual()
    fixedLegConvention = QuantLib.Time.Unadjusted()
    floatingLegConvention = QuantLib.Time.ModifiedFollowing()
    fixedLegDayCounter = QuantLib.Time.EuroThirty360()
    floatingLegFrequency = QuantLib.Time.Semiannual()

    swapType = Payer()
    dummyFixedRate = 0.03
    indexSixMonths = euribor_index(QuantLib.Time.TenorPeriod(Dates.Month(6)), rhTermStructure)

    startDate = QuantLib.Time.advance(Dates.Year(1), cal, settlementDate, floatingLegConvention)
    maturity = QuantLib.Time.advance(Dates.Year(5), cal, startDate, floatingLegConvention)

    fixedSchedule = QuantLib.Time.Schedule(startDate, maturity, QuantLib.Time.TenorPeriod(fixedLegFrequency),
                    fixedLegConvention, fixedLegConvention, QuantLib.Time.DateGenerationForwards(), false, cal)
    floatSchedule = QuantLib.Time.Schedule(startDate, maturity, QuantLib.Time.TenorPeriod(floatingLegFrequency),
                    floatingLegConvention, floatingLegConvention, QuantLib.Time.DateGenerationForwards(), false, cal)

    swap = VanillaSwap(swapType, 1000.0, fixedSchedule, dummyFixedRate, fixedLegDayCounter,
          indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))

    fixedATMRate = fair_rate(swap)

    atmSwap = VanillaSwap(swapType, 1000.0, fixedSchedule, fixedATMRate, fixedLegDayCounter, indexSixMonths,
              0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))

Construct our model

.. code-block:: julia

    modelG2 = G2(rhTermStructure)

Build our calibration helpers

.. code-block:: julia

    numRows = 5
    numCols = 5

    times = zeros(0)
    swaptions = Vector{SwaptionHelper}(numRows)

    for i = 1:numRows
      j = numCols - (i - 1)
      k = (i - 1) * numCols + j

      sh = SwaptionHelper(swaptionMats[i], swaptionLengths[j], Quote(swaptionVols[k]), indexSixMonths,
            indexSixMonths.tenor, indexSixMonths.dc, indexSixMonths.dc, rhTermStructure,
            G2SwaptionEngine(modelG2, 6.0, 16))

      times = add_times_to!(sh, times)
      swaptions[i] = sh
    end

    tg = QuantLib.Time.TimeGrid(times, 30)

Calibrate our model

.. code-block:: julia

    om = QuantLib.Math.LevenbergMarquardt()
    calibrate!(modelG2, swaptions, om, QuantLib.Math.EndCriteria(400, 100, 1.0e-8, 1.0e-8, 1.0e-8))

    for i=1:numRows
      j = numCols - (i - 1)
      k = (i - 1) * numCols + j

      npv = model_value!(swaptions[i])
      implied = implied_volatility!(swaptions[i], npv, 1e-4, 1000, 0.05, 0.50)
      diff = implied - swaptionVols[k]

      println(@sprintf("%i x %i: model %.5f%%, market: %.5f%% (%.5f%%)", i, Int(swaptionLengths[j]), implied * 100,
              swaptionVols[k] * 100, diff * 100))
    end

    println("calibrated to: ")
    println(@sprintf("a = %.6f, sigma = %.6f", get_params(modelG2)[1], get_params(modelG2)[2]))
    println(@sprintf("b = %.6f, eta = %.6f", get_params(modelG2)[3], get_params(modelG2)[4]))
    println(@sprintf("rho = %.6f", get_params(modelG2)[5]))

Build a Bermudan swaption for pricing

.. code-block:: julia

    swapLeg = swap.legs[1] # Fixed Leg

    bermudanDates = Vector{Date}(length(swapLeg.coupons))
    for i=1:length(swapLeg.coupons)
    bermudanDates[i]  = accrual_start_date(swapLeg.coupons[i])
    end

    bermudanExercise = BermudanExercise(bermudanDates)

    bermudanSwaption = Swaption(atmSwap, bermudanExercise)

Use a tree swaption engine to price the swaption with our G2 model

.. code-block:: julia

    bermudanSwaption = update_pricing_engine(bermudanSwaption, TreeSwaptionEngine(modelG2, 50))

    println(@sprintf("G2 (tree):       %.6f", npv(bermudanSwaption)))

Use a finite-differences swaption engine to price the swaption with our G2 model

.. code-block:: julia

    bermudanSwaption = update_pricing_engine(bermudanSwaption, FdG2SwaptionEngine(modelG2))

    println(@sprintf("G2 (fdm):       %.6f", npv(bermudanSwaption)))
