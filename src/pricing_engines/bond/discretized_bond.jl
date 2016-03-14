type DiscretizedDiscountBond <: DiscretizedAsset
  common::DiscretizedAssetCommon
end

DiscretizedDiscountBond() = DiscretizedDiscountBond(DiscretizedAssetCommon())

function reset!(dBond::DiscretizedDiscountBond, sz::Int)
  dBond.common.values = ones(sz)

  return dBond
end

type DiscretizedCallableFixedRateBond <: DiscretizedAsset
  redemptionTime::Float64
  couponTimes::Vector{Float64}
  callabilityTimes::Vector{Float64}
  args::CallableBondArgs
  common::DiscretizedAssetCommon
end

function DiscretizedCallableFixedRateBond(bond::CallableFixedRateBond, referenceDate::Date, dc::DayCount)
  args = CallableBondArgs(bond)
  redemptionTime = year_fraction(dc, referenceDate, args.redemptionDate)
  couponTimes = Vector{Float64}(length(args.couponDates))

  for i in eachindex(couponTimes)
    couponTimes[i] = year_fraction(dc, referenceDate, args.couponDates[i])
  end

  callabilityTimes = Vector{Float64}(length(args.callabilityDates))

  for i in eachindex(callabilityTimes)
    callabilityTimes[i] = year_fraction(dc, referenceDate, args.callabilityDates[i])
  end

  # similar to tree swaption engine, we collapse similar coupon and exercise dates to avoid mispricing.
  # Can remove if necessary

  for i in eachindex(callabilityTimes), j in eachindex(couponTimes)
    if within_next_week(callabilityTimes[i], couponTimes[j])
      couponTimes[j] = callabilityTimes[i]
    end
  end

  return DiscretizedCallableFixedRateBond(redemptionTime, couponTimes, callabilityTimes, args, DiscretizedAssetCommon())
end

function reset!(dcb::DiscretizedCallableFixedRateBond, sz::Int)
  dcb.common.values = fill(dcb.args.redemption, sz)

  adjust_values!(dcb)
  return dcb
end

function mandatory_times(dcb::DiscretizedCallableFixedRateBond)
  times =  dcb.couponTimes[dcb.couponTimes .>= 0.0]
  if dcb.redemptionTime >= 0.0
    insert!(times, 1, dcb.redemptionTime)
  end
  callTimes = dcb.callabilityTimes[dcb.callabilityTimes .>= 0.0]
  times = vcat(times, callTimes)
  return times
end

function apply_callability!(dcb::DiscretizedCallableFixedRateBond, ::Call, i::Int)
  for j in eachindex(dcb.common.values)
    dcb.common.values[j] = min(dcb.args.callabilityPrices[i], dcb.common.values[j])
  end

  return dcb
end

function apply_callability!(dcb::DiscretizedCallableFixedRateBond, ::Put, i::Int)
  for j in eachindex(dcb.common.values)
    dcb.common.values[j] = max(dcb.common.values[j], dcb.args.callabilityPrices[i])
  end

  return dcb
end

function add_coupon!(dcb::DiscretizedCallableFixedRateBond, i::Int)
  dcb.common.values += dcb.args.couponAmounts[i]

  return dcb
end

apply_callability!(dcb::DiscretizedCallableFixedRateBond, i::Int) = apply_callability!(dcb, dcb.args.putCallSchedule[i].optionType, i)

function post_adjust_values_impl!(dcb::DiscretizedCallableFixedRateBond)
  for i in eachindex(dcb.callabilityTimes)
    t = dcb.callabilityTimes[i]
    if t >= 0.0 && is_on_time(dcb, t)
      apply_callability!(dcb, i)
    end
  end

  for i in eachindex(dcb.couponTimes)
    t = dcb.couponTimes[i]
    if t >= 0.0 && is_on_time(dcb, t)
      add_coupon!(dcb, i)
    end
  end

  return dcb
end
