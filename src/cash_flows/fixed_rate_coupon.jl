using QuantLib.Time

mutable struct FixedRateCoupon{DC <: DayCount, IR <: InterestRate} <: Coupon
  couponMixin::CouponMixin{DC}
  paymentDate::Date
  nominal::Float64
  rate::IR
end

FixedRateCoupon(paymentDate::Date,
                faceAmount::Float64,
                rate::IR,
                accrual_start::Date,
                accrual_end::Date,
                ref_start::Date,
                ref_end::Date,
                dc::DC,
                accrual::Float64) where {DC <: DayCount, IR <: InterestRate} =
                FixedRateCoupon{DC, IR}(CouponMixin{DC}(accrual_start, accrual_end, ref_start, ref_end, dc, accrual), paymentDate, faceAmount, rate)

## COUPON METHODS ##
amount(coup::FixedRateCoupon) =
        coup.nominal * (compound_factor(coup.rate, accrual_start_date(coup), accrual_end_date(coup), ref_period_start(coup), ref_period_end(coup)) - 1)

calc_rate(coup::FixedRateCoupon) = coup.rate.rate

mutable struct FixedRateLeg{FRC <: FixedRateCoupon} <: Leg
  coupons::Vector{FRC}
  redemption::Union{SimpleCashFlow, Nothing}
end

function FixedRateLeg(schedule::Schedule,
                      faceAmount::Float64,
                      rate::Float64,
                      calendar::BusinessCalendar,
                      paymentConvention::BusinessDayConvention,
                      dc::DayCount;
                      add_redemption::Bool = true)
  # n = add_redemption ? length(schedule.dates) : length(schedule.dates) - 1
  # coups = Vector{Union{FixedRateCoupon, SimpleCashFlow}}(n)
  #
  # start_date = schedule.dates[1]
  # end_date = schedule.dates[2]
  # payment_date = adjust(calendar, paymentConvention, end_date)
  # #TODO: setup payment adjustments and the like
  # coups[1] = FixedRateCoupon(payment_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)
  #
  # # build coupons
  # count = 2
  # start_date = end_date
  # end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
  # payment_date = adjust(calendar, paymentConvention, end_date)
  # while start_date < schedule.dates[end]
  #   @inbounds coups[count] = FixedRateCoupon(payment_date, faceAmount, InterestRate(rate, dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)
  #
  #   count += 1
  #   start_date = end_date
  #   end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
  #   payment_date = adjust(calendar, paymentConvention, end_date)
  # end
  #
  # if add_redemption
  #   @inbounds coups[end] = SimpleCashFlow(faceAmount, end_date)
  # end

  ratesLen = length(schedule.dates) - 1

  ratesVec = fill(rate, ratesLen)

  FixedRateLeg(schedule, faceAmount, ratesVec, calendar, paymentConvention, dc; add_redemption = add_redemption)
end

function FixedRateLeg(schedule::Schedule,
                      faceAmount::Float64,
                      rates::Vector{Float64},
                      calendar::BusinessCalendar,
                      paymentConvention::BusinessDayConvention,
                      dc::DC;
                      add_redemption::Bool = false) where {DC <: DayCount}
  n = length(schedule.dates) - 1
  length(rates) == length(schedule.dates) - 1 || error("mismatch in coupon rates")
  
  coup_type = FixedRateCoupon{DC, InterestRate{DC, SimpleCompounding, typeof(schedule.tenor.freq)}}
  coups = Vector{coup_type}(undef, n)

  start_date = schedule.dates[1]
  end_date = schedule.dates[2]
  payment_date = adjust(calendar, paymentConvention, end_date)
  #TODO: setup payment adjustments and the like
  coups[1] = FixedRateCoupon(payment_date, faceAmount, InterestRate(rates[1], dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)

  # build coupons
  count = 2
  start_date = end_date
  end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
  payment_date = adjust(calendar, paymentConvention, end_date)
  while start_date < schedule.dates[end]
    @inbounds coups[count] = FixedRateCoupon(payment_date, faceAmount, InterestRate(rates[count], dc, SimpleCompounding(), schedule.tenor.freq), start_date, end_date, start_date, end_date, dc, -1.0)

    count += 1
    start_date = end_date
    end_date = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    payment_date = adjust(calendar, paymentConvention, end_date)
  end

  if add_redemption
    redempt = SimpleCashFlow(faceAmount, end_date)
  else
    redempt = nothing
  end

  return FixedRateLeg{coup_type}(coups, redempt)
end

# Coupon methods
get_pay_dates(coups::Vector{FixedRateCoupon}) = Date[date(coup) for coup in coups]
get_reset_dates(coups::Vector{FixedRateCoupon}) = Date[accrual_start_date(coup) for coup in coups]

function accrued_amount(coup::FixedRateCoupon, settlement_date::Date)
  if settlement_date <= accrual_start_date(coup) || settlement_date > coup.paymentDate
    return 0.0
  end

  return coup.nominal *
      (compound_factor(coup.rate, accrual_start_date(coup), min(settlement_date, accrual_end_date(coup)), ref_period_start(coup), ref_period_end(coup)) - 1.0)
end
