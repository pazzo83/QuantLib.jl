using QuantLib.Time

type FixedRateCoupon{DC <: DayCount} <: Coupon
  couponMixin::CouponMixin{DC}
  paymentDate::Date
  nominal::Float64
  rate::InterestRate
end

FixedRateCoupon{DC <: DayCount}(paymentDate::Date, faceAmount::Float64, rate::InterestRate, accrual_start::Date, accrual_end::Date, ref_start::Date, ref_end::Date, dc::DC, accrual::Float64) =
                FixedRateCoupon(CouponMixin{DC}(accrual_start, accrual_end, ref_start, ref_end, dc, accrual), paymentDate, faceAmount, rate)

## COUPON METHODS ##
amount(coup::FixedRateCoupon) =
        coup.nominal * (compound_factor(coup.rate, accrual_start_date(coup), accrual_end_date(coup), ref_period_start(coup), ref_period_end(coup)) - 1)

calc_rate(coup::FixedRateCoupon) = coup.rate.rate

type FixedRateLeg <: Leg
  coupons::Vector{Union{FixedRateCoupon, SimpleCashFlow}}
  # redemption::SimpleCashFlow

  function FixedRateLeg{B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount}(schedule::Schedule, faceAmount::Float64, rate::Float64, calendar::B, paymentConvention::C, dc::DC; add_redemption::Bool = true)
    n = add_redemption ? length(schedule.dates) : length(schedule.dates) - 1
    coups = Vector{Union{FixedRateCoupon, SimpleCashFlow}}(n)

    start_date = schedule.dates[1]
    end_date = schedule.dates[2]
    payment_date = adjust(calendar, paymentConvention, end_date)
    #TODO: setup payment adjustments and the like
    coups[1] = FixedRateCoupon(payment_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)

    # build coupons
    count = 2
    start_date = end_date
    end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    payment_date = adjust(calendar, paymentConvention, end_date)
    while start_date < schedule.dates[end]
      @inbounds coups[count] = FixedRateCoupon(payment_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)

      count += 1
      start_date = end_date
      end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
      payment_date = adjust(calendar, paymentConvention, end_date)
    end

    if add_redemption
      @inbounds coups[end] = SimpleCashFlow(faceAmount, end_date)
    end

    new(coups)
  end
end

# Coupon methods
get_pay_dates(coups::Vector{Union{FixedRateCoupon, SimpleCashFlow}}) = Date[date(coup) for coup in filter(check_coupon, coups)]
get_reset_dates(coups::Vector{Union{FixedRateCoupon, SimpleCashFlow}}) = Date[accrual_start_date(coup) for coup in filter(check_coupon, coups)]

function accrued_amount(coup::FixedRateCoupon, settlement_date::Date)
  if settlement_date <= accrual_start_date(coup) || settlement_date > coup.paymentDate
    return 0.0
  end

  return coup.nominal *
      (compound_factor(coup.rate, accrual_start_date(coup), min(settlement_date, accrual_end_date(coup)), ref_period_start(coup), ref_period_end(coup)) - 1.0)
end
