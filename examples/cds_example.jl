using QuantLib

generate_flatforward_ts{C <: QuantLib.Time.BusinessCalendar}(cal::C, settlementDate::Date, flatRate::Quote) =
                        FlatForwardTermStructure(settlementDate, cal, flatRate, QuantLib.Time.Actual365())

function main()
  cal = QuantLib.Time.TargetCalendar()
  todays_date = Date(2007, 5, 15)
  settlementDate = todays_date
  set_eval_date!(settings, todays_date)

  flatRate = Quote(0.01)

  tsCurve = generate_flatforward_ts(cal, settlementDate, flatRate)

  recoveryRate = 0.5
  quoteSpreads = [0.0150, 0.0150, 0.0150, 0.0150]
  tenors = [Dates.Month(3), Dates.Month(6), Dates.Year(1), Dates.Year(2)]

  maturities = [QuantLib.Time.adjust(cal, QuantLib.Following(), todays_date + ten) for ten in tenors]

  insts = SpreadCDSHelper[SpreadCDSHelper(Quote(quoteSpreads[i]), tenors[i], 0, cal, QuantLib.Time.Quarterly(), QuantLib.Time.Following(), QuantLib.Time.DateGenerationTwentieth(),
            QuantLib.Time.Actual365(), recoveryRate, tsCurve) for i in eachindex(tenors)]

  hazardRateStructure = PiecewiseDefaultCurve(todays_date, insts, QuantLib.Time.Actual365(), QuantLib.Math.BackwardFlatInterpolation(), HazardRate(), 1.0e-12)

  hr_curve_data = nodes(hazardRateStructure)

  println(@sprintf("1Y Survival Probability: %.6f %%", survival_probability(hazardRateStructure, todays_date + Dates.Year(1)) * 100.0))
  println(@sprintf("2Y Survival Probability: %.6f %%", survival_probability(hazardRateStructure, todays_date + Dates.Year(2)) * 100.0))
end
