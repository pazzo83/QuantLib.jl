using JQuantLib

const numRows = 5
const numCols = 5
const swaptionVols = [0.1490, 0.1340, 0.1228, 0.1189, 0.1148, 0.1290, 0.1201, 0.1146, 0.1108, 0.1040, 0.1149, 0.1112, 0.1070, 0.1010, 0.0957, 0.1047, 0.1021, 0.0980, 0.0951, 0.1270, 0.1000, 0.0950, 0.0900, 0.1230, 0.1160]
const swaptionLengths = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]

function generate_flatforward_ts{C <: JQuantLib.Time.BusinessCalendar}(cal::C, settlementDate::Date)
  flat_rate = Quote(0.04875825)

  ffts = FlatForwardTermStructure(settlementDate, cal, flat_rate, JQuantLib.Time.Actual365())

  return ffts
end

function calibrate_model{M <: ShortRateModel}(model::M, helpers::Vector{SwaptionHelper})
  om = JQuantLib.Math.LevenbergMarquardt()
  calibrate!(model, helpers, om, JQuantLib.Math.EndCriteria(400, 100, 1.0e-8, 1.0e-8, 1.0e-8))

  for i=1:numRows
    j = numCols - (i - 1)
    k = (i - 1) * numCols + j

    npv = model_value!(helpers[i])
    implied = implied_volatility!(helpers[i], npv, 1e-4, 1000, 0.05, 0.50)
    diff = implied - swaptionVols[k]

    println(@sprintf("%i x %i: model %.5f%%, market: %.5f%% (%.5f%%)", i, Int(swaptionLengths[j]), implied * 100, swaptionVols[k] * 100, diff * 100))
  end
end


function main()
  cal = JQuantLib.Time.TargetCalendar()
  settlementDate = Date(2002, 2, 19)
  todays_date = Date(2002, 2, 15)
  set_eval_date!(settings, todays_date)

  const swaptionMats = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]

  # flat yield term strucutre implying 1x5 swap at 5%
  rhTermStructure = generate_flatforward_ts(cal, settlementDate)

  # Define the ATM/OTM/ITM swaps
  fixedLegFrequency = JQuantLib.Time.Annual()
  fixedLegConvention = JQuantLib.Time.Unadjusted()
  floatingLegConvention = JQuantLib.Time.ModifiedFollowing()
  fixedLegDayCounter = JQuantLib.Time.EuroThirty360()
  floatingLegFrequency = JQuantLib.Time.Semiannual()

  swapType = Payer()
  dummyFixedRate = 0.03
  indexSixMonths = euribor_index(JQuantLib.Time.TenorPeriod(Dates.Month(6)), rhTermStructure)

  startDate = JQuantLib.Time.advance(Dates.Year(1), cal, settlementDate, floatingLegConvention)
  maturity = JQuantLib.Time.advance(Dates.Year(5), cal, startDate, floatingLegConvention)

  fixedSchedule = JQuantLib.Time.Schedule(startDate, maturity, JQuantLib.Time.TenorPeriod(fixedLegFrequency), fixedLegConvention, fixedLegConvention, JQuantLib.Time.DateGenerationForwards(), false, cal)
  floatSchedule = JQuantLib.Time.Schedule(startDate, maturity, JQuantLib.Time.TenorPeriod(floatingLegFrequency), floatingLegConvention, floatingLegConvention, JQuantLib.Time.DateGenerationForwards(), false, cal)

  swap = VanillaSwap(swapType, 1000.0, fixedSchedule, dummyFixedRate, fixedLegDayCounter, indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))

  fixedATMRate = fair_rate(swap)
  fixedOTMRate = fixedATMRate * 1.2
  fixedITMRate = fixedATMRate * 0.8

  atmSwap = VanillaSwap(swapType, 1000.0, fixedSchedule, fixedATMRate, fixedLegDayCounter, indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))
  otmSwap = VanillaSwap(swapType, 1000.0, fixedSchedule, fixedOTMRate, fixedLegDayCounter, indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))
  itmSwap = VanillaSwap(swapType, 1000.0, fixedSchedule, fixedITMRate, fixedLegDayCounter, indexSixMonths, 0.0, floatSchedule, indexSixMonths.dc, DiscountingSwapEngine(rhTermStructure))

  times = zeros(0)
  swaptions = Vector{SwaptionHelper}(numRows)

  # First Model
  modelG2 = G2(rhTermStructure)
  hullWhiteModel = HullWhite(rhTermStructure)
  hullWhiteModel2 = HullWhite(rhTermStructure)
  blackKarasinski = BlackKarasinski(rhTermStructure)

  for i = 1:numRows
    j = numCols - (i - 1)
    k = (i - 1) * numCols + j

    sh = SwaptionHelper(swaptionMats[i], swaptionLengths[j], Quote(swaptionVols[k]), indexSixMonths, indexSixMonths.tenor, indexSixMonths.dc, indexSixMonths.dc,
                        rhTermStructure, G2SwaptionEngine(modelG2, 6.0, 16))

    times = add_times_to!(sh, times)
    swaptions[i] = sh
  end

  tg = JQuantLib.Time.TimeGrid(times, 30)

  # models
  # modelG2 = G2(rhTermStructure)
  # hullWhiteModel = HullWhite(rhTermStructure)
  # hullWhiteModel2 = HullWhite(rhTermStructure)
  # blackKarasinski = BlackKarasinski(rhTermStructure)
  #
  # for swaptionHelper in swaptions
  #   swaptionHelper.pricingEngine = G2SwaptionEngine(modelG2, 6.0, 16)
  # end

  println("G2 (analytic formulae) calibration")
  calibrate_model(modelG2, swaptions)
  println("calibrated to: ")
  println(@sprintf("a = %.6f, sigma = %.6f", get_params(modelG2)[1], get_params(modelG2)[2]))
  println(@sprintf("b = %.6f, eta = %.6f", get_params(modelG2)[3], get_params(modelG2)[4]))
  println(@sprintf("rho = %.6f", get_params(modelG2)[5]))

  # Hull White
  for i in eachindex(swaptions)
    swaptions[i] = update_pricing_engine(swaptions[i], JamshidianSwaptionEngine(hullWhiteModel))
  end

  println("")
  println("Hull-White (analytic formulae) calibration")
  calibrate_model(hullWhiteModel, swaptions)
  println("calibrated to: ")
  println(@sprintf("a = %.6f, sigma = %.6f", get_params(hullWhiteModel)[1], get_params(hullWhiteModel)[2]))

  # Hull White (numeric)
  for i in eachindex(swaptions)
    swaptions[i] = update_pricing_engine(swaptions[i], TreeSwaptionEngine(hullWhiteModel2, tg))
  end

  println("")
  println("Hull-White (numerical) calibration")
  calibrate_model(hullWhiteModel2, swaptions)
  println("calibrated to: ")
  println(@sprintf("a = %.6f, sigma = %.6f", get_params(hullWhiteModel2)[1], get_params(hullWhiteModel2)[2]))

  # Black Karasinski
  for i in eachindex(swaptions)
    swaptions[i] = update_pricing_engine(swaptions[i], TreeSwaptionEngine(blackKarasinski, tg))
  end

  println("")
  println("Black-Karasinski (numerical) calibration")
  calibrate_model(blackKarasinski, swaptions)
  println("calibrated to: ")
  println(@sprintf("a = %.6f, sigma = %.6f", get_params(blackKarasinski)[1], get_params(blackKarasinski)[2]))

  # ATM Bermudan swaption pricing
  println("")
  println(@sprintf("Payer bermudan swaption struck at %.6f %% (ATM)", fixedATMRate * 100.0))

  swapLeg = swap.legs[1] # Fixed Leg

  bermudanDates = Vector{Date}(length(swapLeg.coupons))
  for i=1:length(swapLeg.coupons)
    bermudanDates[i]  = accrual_start_date(swapLeg.coupons[i])
  end

  bermudanExercise = BermudanExercise(bermudanDates)

  bermudanSwaption = Swaption(atmSwap, bermudanExercise)

  bermudanSwaption = update_pricing_engine(bermudanSwaption, TreeSwaptionEngine(modelG2, 50))

  println(@sprintf("G2 (tree):       %.6f", npv(bermudanSwaption)))

  bermudanSwaption = update_pricing_engine(bermudanSwaption, FdG2SwaptionEngine(modelG2))

  println(@sprintf("G2 (fdm):       %.6f", npv(bermudanSwaption)))

  bermudanSwaption = update_pricing_engine(bermudanSwaption, TreeSwaptionEngine(hullWhiteModel, 50))

  println(@sprintf("HW (tree):       %.6f", npv(bermudanSwaption)))

  bermudanSwaption = update_pricing_engine(bermudanSwaption, FdHullWhiteSwaptionEngine(hullWhiteModel))

  println(@sprintf("HW (fdm):       %.6f", npv(bermudanSwaption)))

  bermudanSwaption = update_pricing_engine(bermudanSwaption, TreeSwaptionEngine(hullWhiteModel2, 50))

  println(@sprintf("HW (num, tree):       %.6f", npv(bermudanSwaption)))

  bermudanSwaption = update_pricing_engine(bermudanSwaption, FdHullWhiteSwaptionEngine(hullWhiteModel2))

  println(@sprintf("HW (num, fdm):       %.6f", npv(bermudanSwaption)))

  bermudanSwaption = update_pricing_engine(bermudanSwaption, TreeSwaptionEngine(blackKarasinski, 50))

  println(@sprintf("BK:       %.6f", npv(bermudanSwaption)))

  # # OTM Bermudan swaption pricing
  # println("")
  # println(@sprintf("Payer bermudan swaption struck at %.6f %% (OTM)", fixedOTMRate * 100.0))
  #
  # otmBermudanSwaption = Swaption(otmSwap, bermudanExercise)
  #
  # update_pricing_engine!(otmBermudanSwaption, TreeSwaptionEngine(modelG2, 300))
  #
  # println(@sprintf("G2 (tree):       %.6f", npv(otmBermudanSwaption)))
  #
  # update_pricing_engine!(otmBermudanSwaption, FdG2SwaptionEngine(modelG2))
  #
  # println(@sprintf("G2 (fdm):       %.6f", npv(otmBermudanSwaption)))
  #
  # update_pricing_engine!(otmBermudanSwaption, TreeSwaptionEngine(hullWhiteModel, 50))
  #
  # println(@sprintf("HW (tree):       %.6f", npv(otmBermudanSwaption)))
  #
  # update_pricing_engine!(otmBermudanSwaption, FdHullWhiteSwaptionEngine(hullWhiteModel))
  #
  # println(@sprintf("HW (fdm):       %.6f", npv(otmBermudanSwaption)))
  #
  # update_pricing_engine!(otmBermudanSwaption, TreeSwaptionEngine(hullWhiteModel2, 50))
  #
  # println(@sprintf("HW (num, tree):       %.6f", npv(otmBermudanSwaption)))
  #
  # update_pricing_engine!(otmBermudanSwaption, FdHullWhiteSwaptionEngine(hullWhiteModel2))
  #
  # println(@sprintf("HW (num, fdm):       %.6f", npv(otmBermudanSwaption)))
  #
  # update_pricing_engine!(otmBermudanSwaption, TreeSwaptionEngine(blackKarasinski, 50))
  #
  # println(@sprintf("BK:       %.6f", npv(otmBermudanSwaption)))
  #
  # # ITM Bermudan swaption pricing
  # println("")
  # println(@sprintf("Payer bermudan swaption struck at %.6f %% (OTM)", fixedITMRate * 100.0))
  #
  # itmBermudanSwaption = Swaption(itmSwap, bermudanExercise)
  #
  # update_pricing_engine!(itmBermudanSwaption, TreeSwaptionEngine(modelG2, 50))
  #
  # println(@sprintf("G2 (tree):       %.6f", npv(itmBermudanSwaption)))
  #
  # update_pricing_engine!(itmBermudanSwaption, FdG2SwaptionEngine(modelG2))
  #
  # println(@sprintf("G2 (fdm):       %.6f", npv(itmBermudanSwaption)))
  #
  # update_pricing_engine!(itmBermudanSwaption, TreeSwaptionEngine(hullWhiteModel, 50))
  #
  # println(@sprintf("HW (tree):       %.6f", npv(itmBermudanSwaption)))
  #
  # update_pricing_engine!(itmBermudanSwaption, FdHullWhiteSwaptionEngine(hullWhiteModel))
  #
  # println(@sprintf("HW (fdm):       %.6f", npv(itmBermudanSwaption)))
  #
  # update_pricing_engine!(itmBermudanSwaption, TreeSwaptionEngine(hullWhiteModel2, 50))
  #
  # println(@sprintf("HW (num, tree):       %.6f", npv(itmBermudanSwaption)))
  #
  # update_pricing_engine!(itmBermudanSwaption, FdHullWhiteSwaptionEngine(hullWhiteModel2))
  #
  # println(@sprintf("HW (num, fdm):       %.6f", npv(itmBermudanSwaption)))
  #
  # update_pricing_engine!(itmBermudanSwaption, TreeSwaptionEngine(blackKarasinski, 50))
  #
  # println(@sprintf("BK:       %.6f", npv(itmBermudanSwaption)))
end
