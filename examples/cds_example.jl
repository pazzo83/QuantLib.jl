using JQuantLib

generate_flatforward_ts{C <: JQuantLib.Time.BusinessCalendar}(cal::C, settlementDate::Date, flatRate::Quote) =
                        FlatForwardTermStructure(settlementDate, cal, flatRate, JQuantLib.Time.Actual365())

function main()
  cal = JQuantLib.Time.TargetCalendar()
  todays_date = Date(2007, 5, 15)
  settlementDate = todays_date
  set_eval_date!(settings, todays_date)

  flatRate = Quote(0.01)

  tsCurve = generate_flatforward_ts(cal, settlementDate, flatRate)

  recoveryRate = 0.5
  quoteSpreads = [0.0150, 0.0150, 0.0150, 0.0150]
  tenors = [Dates.Month(3), Dates.Month(6), Dates.Year(1), Dates.Year(2)]

  maturities = [JQuantLib.Time.adjust(cal, JQuantLib.Following(), todays_date + ten) for ten in tenors]

  insts = SpreadCDSHelper[SpreadCDSHelper(Quote(quoteSpreads[i]), tenors[i], 0, cal, JQuantLib.Time.Quarterly(), JQuantLib.Time.Following(), JQuantLib.Time.DateGenerationTwentieth(),
            JQuantLib.Time.Actual365(), recoveryRate, tsCurve) for i in eachindex(tenors)]

  hazardRateStructure = PiecewiseDefaultCurve(todays_date, insts, JQuantLib.Time.Actual365(), JQuantLib.Math.BackwardFlatInterpolation(), HazardRate(), 1.0e-12)

  return nodes(hazardRateStructure)
end
