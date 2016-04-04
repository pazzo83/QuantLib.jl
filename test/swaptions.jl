using Base.Test
using QuantLib

cal = QuantLib.Time.TargetCalendar()
settlementDate = Date(2002, 2, 19)
todays_date = Date(2002, 2, 15)
set_eval_date!(settings, todays_date)

swaptionMats = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]
swaptionVols = [0.1490, 0.1340, 0.1228, 0.1189, 0.1148, 0.1290, 0.1201, 0.1146, 0.1108,
                0.1040, 0.1149, 0.1112, 0.1070, 0.1010, 0.0957, 0.1047, 0.1021, 0.0980, 0.0951,
                0.1270, 0.1000, 0.0950, 0.0900, 0.1230, 0.1160]
swaptionLengths = [Dates.Year(1), Dates.Year(2), Dates.Year(3), Dates.Year(4), Dates.Year(5)]

flat_rate = Quote(0.04875825)
rhTermStructure = FlatForwardTermStructure(settlementDate, cal, flat_rate, QuantLib.Time.Actual365())

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

modelG2 = G2(rhTermStructure)

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

om = QuantLib.Math.LevenbergMarquardt()
calibrate!(modelG2, swaptions, om, QuantLib.Math.EndCriteria(400, 100, 1.0e-8, 1.0e-8, 1.0e-8))

for i=1:numRows
  j = numCols - (i - 1)
  k = (i - 1) * numCols + j

  npv = model_value!(swaptions[i])
  implied = implied_volatility!(swaptions[i], npv, 1e-4, 1000, 0.05, 0.50)
  diff = implied - swaptionVols[k]
end

swapLeg = swap.legs[1] # Fixed Leg

bermudanDates = Vector{Date}(length(swapLeg.coupons))
for i=1:length(swapLeg.coupons)
bermudanDates[i]  = accrual_start_date(swapLeg.coupons[i])
end

bermudanExercise = BermudanExercise(bermudanDates)

bermudanSwaption = Swaption(atmSwap, bermudanExercise)

bermudanSwaption = update_pricing_engine(bermudanSwaption, TreeSwaptionEngine(modelG2, 50))

@test npv(bermudanSwaption) == 14.11010300956427
