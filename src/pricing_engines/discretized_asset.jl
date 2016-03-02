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
  if ~QuantLib.Math.close_enough(dAsset.common.time, dAsset.common.latestPreAdjustment)
    pre_adjust_values_impl!(dAsset)
    dAsset.common.latestPreAdjustment = dAsset.common.time
  end

  return dAsset
end

function post_adjust_values!(dAsset::DiscretizedAsset)
  if ~QuantLib.Math.close_enough(dAsset.common.time, dAsset.common.latestPostAdjustment)
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
  return QuantLib.Math.close_enough(grid.times[findfirst(grid.times .>= t)], dAsset.common.time)
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

function initialize!(asset::DiscretizedAsset, lattice::Lattice, t::Float64)
  set_method!(asset, lattice)
  initialize!(lattice, asset, t)
  return asset
end

function reset!(dBond::DiscretizedDiscountBond, sz::Int)
  dBond.common.values = ones(sz)

  return dBond
end
