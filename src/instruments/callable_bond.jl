type CallableFixedRateBond{DC <: DayCount, P <: PricingEngine, P2 <: PricingEngine} <: AbstractCallableBond
  lazyMixin::LazyMixin
  bondMixin::BondMixin{Int}
  faceAmount::Float64
  schedule::Schedule
  cashflows::FixedRateLeg
  dc::DC
  redemption::Float64
  startDate::Date
  pricingEngine::P
  settlementValue::Float64
  putCallSchedule::CallabilitySchedule
  blackEngine::P2
  blackVolQuote::Quote
end

function CallableFixedRateBond{DC <: DayCount, P <: PricingEngine}(settlementDays::Int, faceAmount::Float64, schedule::Schedule, coupons::Union{Vector{Float64}, Float64},
                              accrualDayCounter::DC, paymentConvention::BusinessDayConvention, redemption::Float64, issueDate::Date,
                              putCallSchedule::CallabilitySchedule, pe::P)
  maturityDate = schedule.dates[end]

  # build bond
  isZeroCouponBond = length(coupons) == 1 && is_close(coupons[1], 0.0)

  if ~isZeroCouponBond
    coups = FixedRateLeg(schedule, faceAmount, coupons, schedule.cal, paymentConvention, accrualDayCounter)
  else
    # build redemption CashFlow
    redemptionDate = adjust(schedule.cal, paymentConvention, maturityDate)
    coups = ZeroCouponLeg(SimpleCashFlow(redemption, redemptionDate))
  end

  blackVolQuote = Quote(0.0)
  blackEngine = BlackCallableFixedRateBondEngine(blackVolQuote)

  return CallableFixedRateBond(LazyMixin(), BondMixin{Int}(settlementDays, issueDate, maturityDate), faceAmount, schedule, coups, accrualDayCounter, redemption,
                              schedule.dates[1], pe, 0.0, putCallSchedule, blackEngine, blackVolQuote)
end
