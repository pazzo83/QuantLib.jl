using Dates

mutable struct DiscretizedSwaption{E <: Exercise, L <: Lattice} <: DiscretizedOption
  underlying::DiscretizedSwap
  exercise::E
  exerciseTimes::Vector{Float64}
  fixedPayDates::Vector{Date}
  fixedResetDates::Vector{Date}
  floatingPayDates::Vector{Date}
  floatingResetDates::Vector{Date}
  lastPayment::Float64
  common::DiscretizedAssetCommon{L}
end

function DiscretizedSwaption(swaption::Swaption, referenceDate::Date, dc::DayCount, lattice::L) where {L <: Lattice}
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

  # for i = 1:n
  #   exerciseTimes[i] = year_fraction(dc, referenceDate, dates[i])
  # end
  map!(x -> year_fraction(dc, referenceDate, x), exerciseTimes, dates)

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
  underlying = DiscretizedSwap(nominal, swapT, referenceDate, dc, fixedPayDates, fixedResetDates, floatingPayDates, floatingResetDates, swaption.swap.args, lattice)
  exercise = swaption.exercise

  DiscretizedSwaption{typeof(exercise), L}(underlying, exercise, exerciseTimes, fixedPayDates, fixedResetDates,
                                          floatingPayDates, floatingResetDates, lastPayment, DiscretizedAssetCommon(lattice))
end

function mandatory_times(discretizedSwaption::DiscretizedSwaption)
  times = mandatory_times(discretizedSwaption.underlying)
  times =  times[times .>= 0.0]
  times = vcat(times, discretizedSwaption.exerciseTimes)
  return times
end

function reset!(dSwaption::DiscretizedSwaption, sz::Int)
  initialize!(dSwaption.underlying, dSwaption.common.method, dSwaption.lastPayment)
  general_reset!(dSwaption, sz)
  return dSwaption
end
