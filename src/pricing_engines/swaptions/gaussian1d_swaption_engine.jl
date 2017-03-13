# Gaussian One-D engines #
abstract GaussianProbabilities
type NoneProbabilities <: GaussianProbabilities end
type NaiveProbabilities <: GaussianProbabilities end
type DigitalProbabilities <: GaussianProbabilities end

type Gaussian1DSwaptionEngine{G <: Gaussian1DModel, Y <: YieldTermStructure, P <: GaussianProbabilities} <: PricingEngine
  model::G
  integrationPoints::Int
  stddevs::Float64
  extrapolatePayoff::Bool
  flatPayoffExtrapolation::Bool
  discountCurve::Y
  probabilities::P

  Gaussian1DSwaptionEngine(model, integrationPoints, stddevs, extrapolatePayoff, flatPayoffExtrapolation, discountCurve, probabilities) =
    new(model, integrationPoints, stddevs, extrapolatePayoff, flatPayoffExtrapolation, discountCurve, probabilities)
end

Gaussian1DSwaptionEngine{G <: Gaussian1DModel, Y <: YieldTermStructure, P <: GaussianProbabilities}(model::G, integrationPoints::Int = 64, stddevs::Float64 = 7.0,
                        extrapolatePayoff::Bool = true, flatPayoffExtrapolation::Bool = false, discountCurve::Y = NullYieldTermStructure(), probabilities::P = NoneProbabilities()) =
                        Gaussian1DSwaptionEngine{G, Y, P}(model, integrationPoints, stddevs, extrapolatePayoff, flatPayoffExtrapolation, discountCurve, probabilities)
# methods #
function _calculate!(pe::Gaussian1DSwaptionEngine, swaption::Swaption)
  isa(swaption.delivery, SettlementPhysical) || error("cash-settled swaptions not yet implemented")
  settlement = reference_date(pe.model.ts)

  if swaption.exercise.dates[end] <= settlement
    swaption.results.value = 0.0 # swaption is expired
    return swaption
  end
  idx = length(swaption.exercise.dates)
  # minIdxAlive = searchsortedlast(swaption.exercise.dates, settlement) + 1
  # minIdxAlive = minIdxAlive == 1 ? 2 : minIdxAlive
  minIdxAlive = upper_bound(swaption.exercise.dates, settlement)

  swap = swaption.swap
  optType = isa(swap.swapT, Payer) ? Call() : Put()
  fixedSchedule = swap.fixedSchedule
  floatingSchedule = swap.floatSchedule

  npv0 = zeros(2 * pe.integrationPoints + 1)
  npv1 = zeros(2 * pe.integrationPoints + 1)
  z = y_grid(pe.model, pe.stddevs, pe.integrationPoints)
  p = zeros(length(z))

  # for probability computation
  npvp0 = Vector{Vector{Float64}}()
  npvp1 = Vector{Vector{Float64}}()
  if !isa(pe.probabilities, NoneProbabilities)
    @inbounds @simd for i = 1:(idx - (minIdxAlive - 1) + 1)
      npvTmp0 = zeros(2 * pe.integrationPoints + 1)
      npvTmp1 = zeros(2 * integrationPoints + 1)
      push!(npvp0, npvTmp0)
      push!(npvp1, npvTmp1)
    end
  end
  # end probability computation

  expiry0 = expiry1 = Date()
  expiry0Time = expiry1Time = -1.0

  while true
    if idx == minIdxAlive - 1
      expiry0 = settlement
    else
      expiry0 = swaption.exercise.dates[idx]
    end

    expiry0Time = max(time_from_reference(pe.model.ts, expiry0), 0.0)
    j1 = upper_bound(fixedSchedule.dates, expiry0 - Dates.Day(1))
    k1 = upper_bound(floatingSchedule.dates, expiry0 - Dates.Day(1))

    # for multi-threading / multi-processing, we need to trigger computation here
    # if expiry0 > settlement
    #   for l = k1:length(get_floating_coupons(swaption.swap))
    #     forward_rate(pe.model, swaption.swap.legs[2].coupons[l].fixingDate, expiry0, 0.0, swaption.swap.iborIndex)
    #     zerobond(pe.model, get_floating_pay_dates(swaption.swap)[l], expiry0, 0.0, pe.discountCurve)
    #   end
    #
    #   for l = j1:length(get_fixed_coupons(swaption.swap))
    #     zerobond(pe.model, get_fixed_pay_dates(swaption.swap)[l], expiry0, 0.0, pe.discountCurve)
    #   end
    #
    #   numeraire(pe.model, expiry0Time, 0.0, pe.discountCurve)
    # end

    @inbounds @simd for k = 1:(expiry0 > settlement ? length(npv0) : 1)
      price = 0.0

      if expiry1Time != -1.0
        yg = y_grid(pe.model, pe.stddevs, pe.integrationPoints, expiry1Time, expiry0Time, expiry0 > settlement ? z[k] : 0.0)
        payoff0 = CubicInterpolation(Spline(), Lagrange(), Lagrange(), z, npv1)
        for i in eachindex(yg)
          p[i] = payoff0(yg[i])
        end

        payoff1 = CubicInterpolation(Spline(), Lagrange(), Lagrange(), z, p)
        for i = 1:length(z) - 1
          price += gaussian_shifted_polynomial_integral(pe.model, 0.0, payoff1.c[i], payoff1.b[i], payoff1.a[i], p[i], z[i], z[i], z[i + 1])
        end

        if pe.extrapolatePayoff
          if pe.flatPayoffExtrapolation
            price += gaussian_shifted_polynomial_integral(pe.model, 0.0, 0.0, 0.0, 0.0, p[length(z) - 1], z[end - 1], z[end], 100.0)
            price += gaussian_shifted_polynomial_integral(pe.model, 0.0, 0.0, 0.0, 0.0, p[1], z[1], -100.0, z[1])
          else
            if isa(optType, Call)
              price += gaussian_shifted_polynomial_integral(pe.model, 0.0, payoff1.c[length(z) - 1], payoff1.b[length(z) - 1], payoff1.a[length(z) - 1], p[length(z) - 1], z[end - 1], z[end], 100.0)
            else
              price += gaussian_shifted_polynomial_integral(pe.model, 0.0, payoff1.c[1], payoff1.b[1], payoff1.a[1], p[1], z[1], -100.0, z[1])
            end
          end
        end
      end
      # println("k: $k, price: $price")
      npv0[k] = price
      # for probability computation
      if !isa(pe.probabilities, NoneProbabilities)
        for m = 1:length(npvp0)
          price = 0.0
          if expiry1Time != -1.0
            yg = y_grid(pe.model, pe.stddevs, pe.integrationPoints, expiry1Time, expiry0Time, expiry0 > settlement ? z[k] : 0.0)
            payoff0 = CubicInterpolation(Spline(), Lagrange(), Lagrange(), z, npvp1[m])
            for i in eachindex(yg)
              p[i] = payoff0(yg[i])
            end

            payoff1 = CubicInterpolation(Spline(), Lagrange(), Lagrange(), z, p)
            for i = 1:length(z) - 1
              price += gaussian_shifted_polynomial_integral(pe.model, 0.0, payoff1.c[i], payoff1.b[i], payoff1.a[i], p[i], z[i], z[i], z[i + 1])
            end

            if pe.extrapolatePayoff
              if pe.flatPayoffExtrapolation
                price += gaussian_shifted_polynomial_integral(pe.model, 0.0, 0.0, 0.0, 0.0, p[length(z) - 1], z[end - 1], z[end], 100.0)
                price += gaussian_shifted_polynomial_integral(pe.model, 0.0, 0.0, 0.0, 0.0, p[1], z[1], -100.0, z[1])
              else
                if isa(optType, Call)
                  price += gaussian_shifted_polynomial_integral(pe.model, 0.0, payoff1.c[length(z) - 1], payoff1.b[length(z) - 1], payoff1.a[length(z) - 1], p[length(z) - 1], z[end - 1], z[end], 100.0)
                else
                  price += gaussian_shifted_polynomial_integral(pe.model, 0.0, payoff1.c[1], payoff1.b[1], payoff1.a[1], p[1], z[1], -100.0, z[1])
                end
              end
            end
          end
          npvp0[m][k] = price
        end
      end
      # end probability computation
      if expiry0 > settlement
        floatingLegNPV = 0.0
        for l = k1:length(get_floating_coupons(swaption.swap))
          floatingLegNPV += swaption.swap.nominal * get_floating_accrual_times(swaption.swap)[l] *
            (get_floating_spreads(swaption.swap)[l] + forward_rate(pe.model, swaption.swap.legs[2].coupons[l].fixingDate, expiry0, z[k], swaption.swap.iborIndex)) *
            zerobond(pe.model, get_floating_pay_dates(swaption.swap)[l], expiry0, z[k], pe.discountCurve)
        end
        # println("floatingLeg: ", floatingLegNPV)
        fixedLegNPV = 0.0
        for l = j1:length(get_fixed_coupons(swaption.swap))
          fixedLegNPV += get_fixed_coupons(swaption.swap)[l] * zerobond(pe.model, get_fixed_pay_dates(swaption.swap)[l], expiry0, z[k], pe.discountCurve)
        end
        # println("fixedLeg: ", fixedLegNPV)
        exerciseValue = (isa(optType, Call) ? 1.0 : -1.0) * (floatingLegNPV - fixedLegNPV) / numeraire(pe.model, expiry0Time, z[k], pe.discountCurve)
        # probability computation
        if !isa(pe.probabilities, NoneProbabilities)
          if idx == length(swaption.exercise.dates) # if true we are at the latest date so we init the no call prob
            npvp0[end][k] = isa(pe.probabilities, NaiveProbabilities) ? 1.0 : 1.0 / (zerobond(pe.model, expiry0Time, 0.0, 0.0, pe.discountCurve) * numeraire(pe.model, expiry0, z[k], pe.discountCurve))
          end

          if exerciseValue >= npv0[k]
            npvp0[idx - (minIdxAlive - 1)][k] = isa(pe.probabilities, NaiveProbabilities) ? 1.0 : (zerobond(pe.model, expiry0Time, 0.0, 0.0, pe.discountCurve) * numeraire(pe.model, expiry0Time, z[k], pe.discountCurve))
            for ii = idx - (minIdxAlive - 1) + 1:length(npvp0)
              npvp0[ii][k] = 0.0
            end
          end
        end
        # end probability computation
        npv0[k] = max(npv0[k], exerciseValue)
      end
    end
    npv1, npv0 = npv0, npv1

    # for probability computation
    if !isa(pe.probabilities, NoneProbabilities)
      for i in eachindex(npvp0)
        npvp1[i], npvp0[i] = npvp0[i], npvp1[i]
      end
    end
    # end probability computation

    expiry1 = expiry0
    expiry1Time = expiry0Time
    idx -= 1
    if idx < minIdxAlive - 1
      break
    end
  end
  # println("NPV1: ", npv1)

  swaption.results.value = npv1[1] * numeraire(pe.model, 0.0, 0.0, pe.discountCurve)

  # for probability computation
  if !isa(pe.probabilities, NoneProbabilities)
    prob = zeros(length(npvp0))
    for i in eachindex(npvp0)
      prob[i] = npvp[i][1] * (isa(pe.probabilities, NaiveProbabilities) ? 1.0 : numeraire(pe.model, 0.0, 0.0, pe.discountCurve))
    end
    swaption.results.additionalResults["probabilities"] = prob
  end

  return swaption
end
