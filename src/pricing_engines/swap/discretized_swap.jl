type DiscretizedSwap{ST <: SwapType} <: DiscretizedAsset
  nominal::Float64
  swapT::ST
  fixedResetTimes::Vector{Float64}
  fixedPayTimes::Vector{Float64}
  floatingResetTimes::Vector{Float64}
  floatingPayTimes::Vector{Float64}
  args::VanillaSwapArgs
  common::DiscretizedAssetCommon
end

function DiscretizedSwap{DC <: DayCount, ST <: SwapType}(nominal::Float64, swapT::ST, referenceDate::Date, dc::DC, fixedPayDates::Vector{Date}, fixedResetDates::Vector{Date},
                        floatingPayDates::Vector{Date}, floatingResetDates::Vector{Date}, args::VanillaSwapArgs)
  fixed_n = length(fixedPayDates)
  float_n = length(floatingPayDates)

  fixedResetTimes = zeros(fixed_n)
  fixedPayTimes = zeros(fixed_n)
  floatingResetTimes = zeros(float_n)
  floatingPayTimes = zeros(float_n)

  for i = 1:fixed_n
    fixedResetTimes[i] = year_fraction(dc, referenceDate, fixedResetDates[i])
    fixedPayTimes[i] = year_fraction(dc, referenceDate, fixedPayDates[i])
  end

  for i = 1:float_n
    floatingResetTimes[i] = year_fraction(dc, referenceDate, floatingResetDates[i])
    floatingPayTimes[i] = year_fraction(dc, referenceDate, floatingPayDates[i])
  end

  DiscretizedSwap(nominal, swapT, fixedResetTimes, fixedPayTimes, floatingResetTimes, floatingPayTimes, args, DiscretizedAssetCommon())
end

# methods #
function mandatory_times(discretizedSwap::DiscretizedSwap)
  # get times
  times = vcat(discretizedSwap.fixedResetTimes[discretizedSwap.fixedResetTimes .>= 0.0], discretizedSwap.fixedPayTimes[discretizedSwap.fixedPayTimes .>= 0.0],
               discretizedSwap.floatingResetTimes[discretizedSwap.floatingResetTimes .>= 0.0], discretizedSwap.floatingPayTimes[discretizedSwap.floatingPayTimes .>= 0.0])

  return times
end

function reset!(dSwap::DiscretizedSwap, sz::Int)
  dSwap.common.values = zeros(sz)
  adjust_values!(dSwap)

  return dSwap
end

function pre_adjust_values_impl!(dSwap::DiscretizedSwap)
  # Floating payments
  @simd for i = 1:length(dSwap.floatingResetTimes)
    @inbounds t = dSwap.floatingResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t)
      bond = DiscretizedDiscountBond()
      initialize!(bond, dSwap.common.method, dSwap.floatingPayTimes[i])
      rollback!(bond, dSwap.common.time)

      nominal = dSwap.nominal
      T = dSwap.args.floatingAccrualTimes[i]
      spread = dSwap.args.floatingSpreads[i]

      accruedSpread = nominal * T * spread
      for j = 1:length(dSwap.common.values)
        @inbounds coup = nominal * (1.0 - bond.common.values[j]) + accruedSpread * bond.common.values[j]

        if isa(dSwap.swapT, Payer)
          @inbounds dSwap.common.values[j] += coup
        else
          @inbounds dSwap.common.values[j] -= coup
        end
      end
    end
  end

  # Fixed Payments
  @simd for i = 1:length(dSwap.fixedResetTimes)
    t = dSwap.fixedResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t)
      bond = DiscretizedDiscountBond()
      initialize!(bond, dSwap.common.method, dSwap.fixedPayTimes[i])
      rollback!(bond, dSwap.common.time)

      @inbounds fixedCoup = dSwap.args.fixedCoupons[i]

      for j = 1:length(dSwap.common.values)
        @inbounds coup = fixedCoup * bond.common.values[j]
        if isa(dSwap.swapT, Payer)
          @inbounds dSwap.common.values[j] -= coup
        else
          @inbounds dSwap.common.values[j] += coup
        end
      end
    end
  end

  return dSwap
end

function post_adjust_values_impl!(dSwap::DiscretizedSwap)
  # fixed coupons whose reset time is in the past won't be managed in pre_adjust_values
  @simd for i = 1:length(dSwap.fixedPayTimes)
    @inbounds t = dSwap.fixedPayTimes[i]
    @inbounds _reset = dSwap.fixedResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t) && _reset < 0.0
      fixedCoup = dSwap.args.fixedCoupons[i]
      if isa(dSwap.swapT, Payer)
        dSwap.common.values -= fixedCoup
      else
        dSwap.common.values += fixedCoup
      end
    end
  end

  # the same applies to floating payments whose rate is already fixed
  @simd for i = 1:length(dSwap.floatingPayTimes)
    @inbounds t = dSwap.floatingPayTimes[i]
    @inbounds _reset = dSwap.floatingResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t) && _reset < 0.0
      @inbounds currentFloatingCoup = dSwap.args.floatingCoupons[i]

      if isa(dSwap.swapT, Payer)
        dSwap.common.values += currentFloatingCoup
      else
        dSwap.common.values -= currentFloatingCoup
      end
    end
  end

  return dSwap
end
