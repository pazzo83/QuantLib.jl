type DiscretizedSwap{ST <: SwapType, L <: Lattice} <: DiscretizedAsset
  nominal::Float64
  swapT::ST
  fixedResetTimes::Vector{Float64}
  fixedPayTimes::Vector{Float64}
  floatingResetTimes::Vector{Float64}
  floatingPayTimes::Vector{Float64}
  args::VanillaSwapArgs
  common::DiscretizedAssetCommon{L}
end

function DiscretizedSwap{DC <: DayCount, ST <: SwapType, L <: Lattice}(nominal::Float64, swapT::ST, referenceDate::Date, dc::DC, fixedPayDates::Vector{Date}, fixedResetDates::Vector{Date},
                        floatingPayDates::Vector{Date}, floatingResetDates::Vector{Date}, args::VanillaSwapArgs, lattice::L)
  fixed_n = length(fixedPayDates)
  float_n = length(floatingPayDates)

  fixedResetTimes = zeros(fixed_n)
  fixedPayTimes = zeros(fixed_n)
  floatingResetTimes = zeros(float_n)
  floatingPayTimes = zeros(float_n)

  @simd for i = 1:fixed_n
    @inbounds fixedResetTimes[i] = year_fraction(dc, referenceDate, fixedResetDates[i])
    @inbounds fixedPayTimes[i] = year_fraction(dc, referenceDate, fixedPayDates[i])
  end

  @simd for i = 1:float_n
    @inbounds floatingResetTimes[i] = year_fraction(dc, referenceDate, floatingResetDates[i])
    @inbounds floatingPayTimes[i] = year_fraction(dc, referenceDate, floatingPayDates[i])
  end

  DiscretizedSwap{ST, L}(nominal, swapT, fixedResetTimes, fixedPayTimes, floatingResetTimes, floatingPayTimes, args, DiscretizedAssetCommon(lattice))
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

_pre_adjust_calc_floating(::Payer, val::Float64, addIt::Float64) = val + addIt
_pre_adjust_calc_floating(::Receiver, val::Float64, addIt::Float64) = val - addIt
_pre_adjust_calc_fixed(::Payer, val::Float64, addIt::Float64) = val - addIt
_pre_adjust_calc_fixed(::Receiver, val::Float64, addIt::Float64) = val + addIt

function pre_adjust_values_impl!(dSwap::DiscretizedSwap)
  # Floating payments
  for i in eachindex(dSwap.floatingResetTimes)
    @inbounds t = dSwap.floatingResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t)
      bond = DiscretizedDiscountBond(dSwap.common.method)
      initialize!(bond, dSwap.common.method, dSwap.floatingPayTimes[i])
      rollback!(bond, dSwap.common.time)

      nominal = dSwap.nominal
      T = dSwap.args.floatingAccrualTimes[i]
      spread = dSwap.args.floatingSpreads[i]

      accruedSpread = nominal * T * spread
      @simd for j in eachindex(dSwap.common.values)
        @inbounds coup = nominal * (1.0 - bond.common.values[j]) + accruedSpread * bond.common.values[j]

        # if isa(dSwap.swapT, Payer)
        #   @inbounds dSwap.common.values[j] += coup
        # else
        #   @inbounds dSwap.common.values[j] -= coup
        @inbounds dSwap.common.values[j] = _pre_adjust_calc_floating(dSwap.swapT, dSwap.common.values[j], coup)
      end
    end
  end

  # Fixed Payments
  for i in eachindex(dSwap.fixedResetTimes)
    t = dSwap.fixedResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t)
      bond = DiscretizedDiscountBond(dSwap.common.method)
      initialize!(bond, dSwap.common.method, dSwap.fixedPayTimes[i])
      rollback!(bond, dSwap.common.time)

      @inbounds fixedCoup = dSwap.args.fixedCoupons[i]

      @simd for j in eachindex(dSwap.common.values)
        # @inbounds coup = fixedCoup * bond.common.values[j]
        # if isa(dSwap.swapT, Payer)
        #   @inbounds dSwap.common.values[j] -= coup
        # else
        #   @inbounds dSwap.common.values[j] += coup
        # end
        @inbounds dSwap.common.values[j] = _pre_adjust_calc_fixed(dSwap.swapT, dSwap.common.values[j], fixedCoup * bond.common.values[j])
      end
    end
  end

  return dSwap
end

_post_adjust_calc_floating(::Payer, val::Vector{Float64}, addIt::Float64) = val + addIt
_post_adjust_calc_floating(::Receiver, val::Vector{Float64}, addIt::Float64) = val - addIt
_post_adjust_calc_fixed(::Payer, val::Vector{Float64}, addIt::Float64) = val - addIt
_post_adjust_calc_fixed(::Receiver, val::Vector{Float64}, addIt::Float64) = val + addIt

function post_adjust_values_impl!(dSwap::DiscretizedSwap)
  # fixed coupons whose reset time is in the past won't be managed in pre_adjust_values
  @simd for i in eachindex(dSwap.fixedPayTimes)
    @inbounds t = dSwap.fixedPayTimes[i]
    @inbounds _reset = dSwap.fixedResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t) && _reset < 0.0
      fixedCoup = dSwap.args.fixedCoupons[i]
      # if isa(dSwap.swapT, Payer)
      #   dSwap.common.values -= fixedCoup
      # else
      #   dSwap.common.values += fixedCoup
      # end
      @inbounds dSwap.common.values = _post_adjust_calc_fixed(dSwap.swapT, dSwap.common.values, dSwap.args.fixedCoupons[i])
    end
  end

  # the same applies to floating payments whose rate is already fixed
  @simd for i in eachindex(dSwap.floatingPayTimes)
    @inbounds t = dSwap.floatingPayTimes[i]
    @inbounds _reset = dSwap.floatingResetTimes[i]
    if t >= 0.0 && is_on_time(dSwap, t) && _reset < 0.0
      # @inbounds currentFloatingCoup = dSwap.args.floatingCoupons[i]

      # if isa(dSwap.swapT, Payer)
      #   dSwap.common.values += currentFloatingCoup
      # else
      #   dSwap.common.values -= currentFloatingCoup
      # end
      @inbounds dSwap.common.values = _post_adjust_calc_floating(dSwap.swapT, dSwap.common.values, dSwap.args.floatingCoupons[i])
    end
  end

  return dSwap
end
