using QuantLib
using Dates

function make_swap(startDate::Date, maturity::Dates.Period, nominal::Float64, fixedRate::Float64, iborIndex::IborIndex, yts::Y, typ::SwapType = Payer()) where {Y <: YieldTermStructure}
  endDate = QuantLib.Time.advance(maturity, QuantLib.Time.TargetCalendar(), startDate)
  fixedLegTenor = Dates.Year(1)
  fixedLegBDC = QuantLib.Time.ModifiedFollowing()
  fixedLegDC = QuantLib.Time.BondThirty360()
  spread = 0.0
  fixedSchedule = QuantLib.Time.Schedule(startDate, endDate, QuantLib.Time.TenorPeriod(fixedLegTenor), fixedLegBDC, fixedLegBDC,
                  QuantLib.Time.DateGenerationBackwards(), false, iborIndex.fixingCalendar)

  floatSchedule = QuantLib.Time.Schedule(startDate, endDate, iborIndex.tenor, iborIndex.convention, iborIndex.convention,
                  QuantLib.Time.DateGenerationBackwards(), false, iborIndex.fixingCalendar)

  vanillaSwap = VanillaSwap(typ, nominal, fixedSchedule, fixedRate, fixedLegDC, iborIndex, spread, floatSchedule, iborIndex.dc, DiscountingSwapEngine{Y}(yts))

  return vanillaSwap, Date[fixing_date(iborIndex, x) for x in floatSchedule.dates[1:end-1]]
end

function main()
  tdy = Date(2015, 4, 7)
  set_eval_date!(settings, tdy)

  rate = 0.03

  dc = QuantLib.Time.Actual365()
  yts = FlatForwardTermStructure(tdy, rate, dc)

  euribor6m = euribor_index(QuantLib.Time.TenorPeriod(Dates.Month(6)), yts)

  portfolio = Tuple{VanillaSwap, Vector{Date}}[make_swap(tdy + Dates.Day(2), Dates.Year(5), 1e6, 0.03, euribor6m, yts),
              make_swap(tdy + Dates.Day(2), Dates.Year(4), 5e5, 0.03, euribor6m, yts, Receiver())]

  for (deal, fixingDates) in portfolio
    npv(deal) # force calculation
  end

  # Stochastic Process
  volas = [0.0075, 0.0075]
  meanRev = 0.02

  model = GSR(yts, [tdy + Dates.Day(100)], volas, meanRev, 16.0)

  process = model.stateProcess

  date_grid = Date[tdy + Dates.Month(i) for i = 0:12*6-1]
  for (_, d) in portfolio
    append!(date_grid, d)
  end
  date_grid = unique(sort(date_grid))
  time_grid = Float64[QuantLib.Time.year_fraction(QuantLib.Time.ISDAActualActual(), tdy, x) for x in date_grid]
  dt = time_grid[2:end] - time_grid[1:end-1]

  seed = 1
  generator = QuantLib.Math.InverseCumulativeRSG(seed, length(time_grid) - 1)

  # Generate N paths
  N = 5000
  x = zeros(length(time_grid), N)
  y = zeros(length(time_grid), N)
  pillars = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
  zero_bonds = zeros(12, length(time_grid), N)

  for j = 1:12
    zero_bonds[j, 1, :] = fill(QuantLib.zerobond(model, pillars[j], 0.0, 0.0, model.ts), N)
  end

  for n = 1:N
    dWs = QuantLib.Math.next_sequence!(generator)[1]
    for i = 2:length(time_grid)
      t0 = time_grid[i-1]
      t1 = time_grid[i]
      x[i, n] = QuantLib.expectation(process, t0, x[i-1, n], dt[i-1]) + dWs[i-1] * QuantLib.std_deviation(process, t0, x[i-1, n], dt[i-1])
      y[i, n] = (x[i, n] - QuantLib.expectation(process, 0.0, 0.0, t1)) / QuantLib.std_deviation(process, 0.0, 0.0, t1)
      for j = 1:12
        zero_bonds[j, i, n] = QuantLib.zerobond(model, t1 + pillars[j], t1, y[i, n], model.ts)
      end
    end
  end

  npv_cube = zeros(length(portfolio), length(date_grid), N)
  # npv_cube = zeros(length(date_grid), length(portfolio), N)
  discount_factors = Float64[discount(yts, time_grid[i]) for i in eachindex(time_grid)]
  # discount_factors = Float64[discount(yts, time_grid[i]) for i in eachindex(time_grid)]

  # Clone swaptions for updating pricing
  tempPort = Tuple{VanillaSwap, Vector{Date}}[]
  defaultDiscCurve = InterpolatedDiscountCurve(Date[tdy for i = 1:12], ones(12), QuantLib.Time.Actual365(), QuantLib.Math.LogLinear())
  for (deal, dates) in portfolio
    push!(tempPort, (QuantLib.clone(deal, DiscountingSwapEngine{typeof(defaultDiscCurve)}(defaultDiscCurve), defaultDiscCurve), dates))
  end
  portfolio = tempPort
  for p = 1:N
    for t in eachindex(date_grid)
      d = date_grid[t]
      set_eval_date!(settings, d)
      ycDates = [d, d + Dates.Month(6)]
      append!(ycDates, [d + Dates.Year(i) for i = 1:10])
      yc = InterpolatedDiscountCurve(ycDates, zero_bonds[:, t, p], QuantLib.Time.Actual365(), QuantLib.Math.LogLinear())
      for deal in portfolio
        update_all_ts!(deal[1], yc)
      end
      if is_valid_fixing_date(euribor6m, d)
        fix = fixing(euribor6m, yc, d)
        add_fixing!(euribor6m, d, fix)
      end

      for i in eachindex(portfolio)
        npv_cube[i, t, p] = npv(portfolio[i][1])
        # npv_cube[t, i, p] = npv(portfolio[i][1])
      end
    end
    empty!(euribor6m.pastFixings)
  end
  set_eval_date!(settings, tdy)
  discountedCube = zeros(size(npv_cube))
  for i = 1:size(discountedCube)[1]
  # for i = 1:size(discountedCube)[2]
    discountedCube[i, :, :] = npv_cube[i, :, :] .* discount_factors
    # discountedCube[:, i, :] = npv_cube[:, i, :] .* discount_factors
  end

  # calculate the discounted rows
  portfolio_npv = sum(npv_cube; dims=1)
  # portfolio_npv = sum(npv_cube, 2)
  discounted_npv = sum(discountedCube; dims=1)
  # discounted_npv = sum(discountedCube, 2)

  # calculate the exposure and discounted exposure
  E = deepcopy(portfolio_npv)
  dE = deepcopy(discounted_npv)

  E[E .< 0.0] .= 0.0
  dE[dE .< 0.0] .= 0.0

  # expected exposure
  EE = sum(E; dims=3) / N
  dEE = sum(dE; dims=3) / N

  PFE_curve = sort(E; dims=3)[:, :, round(Int, 0.95 * N)]
  MPFE = maximum(PFE_curve)

  # alternative pfc 95% quantile of the maxima of each exposure path
  PFE = sort(maximum(E; dims=2), dims=3)[round(Int, 0.95 * N)]
  # PFE = sort(maximum(E, 1), 3)[round(Int, 0.95 * N)]

  # setup default curve
  pdDates = Date[tdy + Dates.Year(i) for i = 0:10]
  hzrates = Float64[0.02 * i for i = 0:10]
  pdCurve = InterpolatedHazardRateCurve(pdDates, hzrates, QuantLib.Time.Actual365(), QuantLib.Math.BackwardFlatInterpolation())

  # calculate defalut probs on grid *times*
  times = collect(range(0, stop=30, length=100))
  dp = zeros(length(times))
  sp = zeros(length(times))
  dd = zeros(length(times))
  hr = zeros(length(times))

  for i in eachindex(times)
    dp[i] = default_probability(pdCurve, times[i])
    sp[i] = survival_probability(pdCurve, times[i])
    dd[i] = default_density(pdCurve, times[i])
    hr[i] = hazard_rate(pdCurve, times[i])
  end

  # Calculation of default probs
  dPD = zeros(length(time_grid) - 1)
  for i = 2:length(time_grid)
    dPD[i-1] = default_probability(pdCurve, time_grid[i-1], time_grid[i])
  end

  # Calculation of CVA
  recovery = 0.4
  CVA = (1.0 - recovery) * sum(dEE[:, 2:end, 1]' .* dPD)
  # CVA = (1.0 - recovery) * sum(dEE[2:end, :, 1] .* dPD)
end
