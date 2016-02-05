## Discretized Asset types ond methods ##

type DiscretizedAssetCommon{L <: Lattice}
  time::Float64
  values::Vector{Float64}
  latestPreAdjustment::Float64
  latestPostAdjustment::Float64
  method::L

  call(::Type{DiscretizedAssetCommon}, t::Float64, v::Vector{Float64}, lpa::Float64, lpoa::Float64) =
      new{Lattice}(t, v, lpa, lpoa)
end

DiscretizedAssetCommon() = DiscretizedAssetCommon(0.0, zeros(0), eps(), eps())

set_time!(a::DiscretizedAsset, t::Float64) = a.common.time = t
set_method!(a::DiscretizedAsset, method::Lattice) = a.common.method = method

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

function DiscretizedSwap{DC <: DayCount, ST <: SwapType}(nominal::Float64, swapT::ST, referenceDate::Date, dc::DC, fixedPayDates::Vector{Date}, fixedResetDates::Vector{Date}, floatingPayDates::Vector{Date}, floatingResetDates::Vector{Date}, args::VanillaSwapArgs)
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

type DiscretizedSwaption{E <: Exercise} <: DiscretizedOption
  underlying::DiscretizedSwap
  exercise::E
  exerciseTimes::Vector{Float64}
  fixedPayDates::Vector{Date}
  fixedResetDates::Vector{Date}
  floatingPayDates::Vector{Date}
  floatingResetDates::Vector{Date}
  lastPayment::Float64
  common::DiscretizedAssetCommon
end

function DiscretizedSwaption{DC <: DayCount}(swaption::Swaption, referenceDate::Date, dc::DC)
  dates = copy(swaption.exercise.dates)
  fixed_coups = swaption.swap.legs[1].coupons
  floating_coups = swaption.swap.legs[2].coupons
  nominal = swaption.swap.nominal
  swapT = swaption.swap.swapT
  n = length(dates)

  exerciseTimes = zeros(n)
  # fixedPayDates = get_pay_dates(fixed_coups)
  # fixedResetDates = get_reset_dates(fixed_coups)
  # floatingPayDates = get_pay_dates(floating_coups)
  # floatingResetDates = get_reset_dates(floating_coups)
  fixedPayDates = copy(swaption.swap.args.fixedPayDates)
  fixedResetDates = copy(swaption.swap.args.fixedResetDates)
  floatingPayDates = copy(swaption.swap.args.floatingPayDates)
  floatingResetDates = copy(swaption.swap.args.floatingResetDates)

  for i = 1:n
    exerciseTimes[i] = year_fraction(dc, referenceDate, dates[i])
  end

  # Date adjustments can get time vectors out of sync
  # Here we try and collapse similar dates which could cause a mispricing
  @simd for i = 1:n
    @inbounds exerciseDate = dates[i]

    for j = 1:length(fixed_coups)
      @inbounds if within_next_week(exerciseDate, fixedPayDates[j]) && fixedResetDates[j] < referenceDate
        @inbounds fixedPayDates[j] = exerciseDate
      end

      @inbounds if within_previous_week(exerciseDate, fixedResetDates[j])
        @inbounds fixedResetDates[j] = exerciseDate
      end
    end

    for j = 1:length(floating_coups)
      @inbounds if within_previous_week(exerciseDate, floatingResetDates[j])
        @inbounds floatingResetDates[j] = exerciseDate
      end
    end
  end

  lastFixedPayment = year_fraction(dc, referenceDate, fixedPayDates[end])
  lastFloatingPayment = year_fraction(dc, referenceDate, floatingPayDates[end])

  lastPayment = max(lastFixedPayment, lastFloatingPayment)
  underlying = DiscretizedSwap(nominal, swapT, referenceDate, dc, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates, swaption.swap.args)
  exercise = swaption.exercise

  DiscretizedSwaption(underlying, exercise, exerciseTimes, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates, lastPayment, DiscretizedAssetCommon())
end

type DiscretizedDiscountBond <: DiscretizedAsset
  common::DiscretizedAssetCommon
end

DiscretizedDiscountBond() = DiscretizedDiscountBond(DiscretizedAssetCommon())

function adjust_values!(dAsset::DiscretizedAsset)
  pre_adjust_values!(dAsset)
  post_adjust_values!(dAsset)

  return dAsset
end

function pre_adjust_values!(dAsset::DiscretizedAsset)
  if ~JQuantLib.Math.close_enough(dAsset.common.time, dAsset.common.latestPreAdjustment)
    pre_adjust_values_impl!(dAsset)
    dAsset.common.latestPreAdjustment = dAsset.common.time
  end

  return dAsset
end

function post_adjust_values!(dAsset::DiscretizedAsset)
  if ~JQuantLib.Math.close_enough(dAsset.common.time, dAsset.common.latestPostAdjustment)
    post_adjust_values_impl!(dAsset)
    dAsset.common.latestPostAdjustment = dAsset.common.time
  end

  return dAsset
end

function partial_rollback!(dAsset::DiscretizedAsset, t::Float64)
  partial_rollback!(dAsset.common.method, dAsset, t)
  return dAsset
end

function is_on_time(dAsset::DiscretizedAsset, t::Float64)
  grid = dAsset.common.method.tg
  return JQuantLib.Math.close_enough(grid.times[findfirst(grid.times .>= t)], dAsset.common.time)
end

function rollback!(dAsset::DiscretizedAsset, t::Float64)
  rollback!(dAsset.common.method, dAsset, t)

  return dAsset
end

pre_adjust_values_impl!(dAsset::DiscretizedAsset) = dAsset # do nothing
post_adjust_values_impl!(dAsset::DiscretizedAsset) = dAsset # do nothing

function general_reset!(dOption::DiscretizedOption, sz::Int)
  # check methods in main and underlying
  dOption.common.values = zeros(sz)
  adjust_values!(dOption)
  return dOption
end

present_value(dAsset::DiscretizedAsset) = present_value(dAsset.common.method, dAsset)


function post_adjust_values_impl!(dOption::DiscretizedOption)
  partial_rollback!(dOption.underlying, dOption.common.time)
  pre_adjust_values!(dOption.underlying)

  if isa(dOption.exercise, AmericanExercise)
    if dOption.common.time >= dOption.exerciseTimes[1] && dOption.common.time <= dOption.exerciseTimes[1]
      apply_exercise_condition!(dOption)
    end
  else
    @simd for i = 1:length(dOption.exerciseTimes)
      @inbounds t = dOption.exerciseTimes[i]
      if t >= 0.0 && is_on_time(dOption, t)
        apply_exercise_condition!(dOption)
      end
    end
  end

  post_adjust_values!(dOption.underlying)

  return dOption
end

function apply_exercise_condition!(dOption::DiscretizedOption)
  @simd for i = 1:length(dOption.common.values)
    @inbounds dOption.common.values[i] = max(dOption.underlying.common.values[i], dOption.common.values[i])
  end

  return dOption
end

function mandatory_times(discretizedSwap::DiscretizedSwap)
  # get times
  times = vcat(discretizedSwap.fixedResetTimes[discretizedSwap.fixedResetTimes .>= 0.0], discretizedSwap.fixedPayTimes[discretizedSwap.fixedPayTimes .>= 0.0],
               discretizedSwap.floatingResetTimes[discretizedSwap.floatingResetTimes .>= 0.0], discretizedSwap.floatingPayTimes[discretizedSwap.floatingPayTimes .>= 0.0])

  return times
end

function mandatory_times(discretizedSwaption::DiscretizedSwaption)
  times = mandatory_times(discretizedSwaption.underlying)
  times =  times[times .>= 0.0]
  times = vcat(times, discretizedSwaption.exerciseTimes)
  return times
end

function initialize!(asset::DiscretizedAsset, lattice::Lattice, t::Float64)
  set_method!(asset, lattice)
  initialize!(lattice, asset, t)
  return asset
end

function reset!(dSwaption::DiscretizedSwaption, sz::Int)
  initialize!(dSwaption.underlying, dSwaption.common.method, dSwaption.lastPayment)
  general_reset!(dSwaption, sz)
  return dSwaption
end

function reset!(dSwap::DiscretizedSwap, sz::Int)
  dSwap.common.values = zeros(sz)
  adjust_values!(dSwap)

  return dSwap
end

function reset!(dBond::DiscretizedDiscountBond, sz::Int)
  dBond.common.values = ones(sz)

  return dBond
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
