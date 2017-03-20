type CallableFixedRateBond{L <: Leg, DC <: DayCount, P <: PricingEngine, P2 <: PricingEngine} <: AbstractCallableBond
  lazyMixin::LazyMixin
  bondMixin::BondMixin
  faceAmount::Float64
  schedule::Schedule
  cashflows::L
  dc::DC
  redemption::Float64
  startDate::Date
  pricingEngine::P
  settlementValue::Float64
  putCallSchedule::CallabilitySchedule
  blackEngine::P2
  blackVolQuote::Quote
end

function CallableFixedRateBond{DC <: DayCount, P <: PricingEngine}(settlementDays::Int, faceAmount::Float64, schedule::Schedule, coupons::Vector{Float64},
                              accrualDayCounter::DC, paymentConvention::BusinessDayConvention, redemption::Float64, issueDate::Date,
                              putCallSchedule::CallabilitySchedule, pe::P)
  maturityDate = schedule.dates[end]

  # build bond
  isZeroCouponBond = length(coupons) == 1 && is_close(coupons[1], 0.0)

  if ~isZeroCouponBond
    couprates = length(coupons) == 1 ? coupons[1] : coupons
    coups = FixedRateLeg(schedule, faceAmount, couprates, schedule.cal, paymentConvention, accrualDayCounter)
  else
    # build redemption CashFlow
    redemptionDate = adjust(schedule.cal, paymentConvention, maturityDate)
    coups = ZeroCouponLeg(SimpleCashFlow(redemption, redemptionDate))
  end

  blackVolQuote = Quote(0.0)
  blackEngine = BlackCallableFixedRateBondEngine(blackVolQuote)

  return CallableFixedRateBond{typeof(coups), DC, P, BlackCallableFixedRateBondEngine}(LazyMixin(), BondMixin(settlementDays, issueDate, maturityDate), faceAmount, schedule, coups, accrualDayCounter, redemption,
                              schedule.dates[1], pe, 0.0, putCallSchedule, blackEngine, blackVolQuote)
end

function accrued(bond::CallableFixedRateBond, d::Date)
  if d == Date()
    d = settlement_date(bond)
  end

  for i in eachindex(bond.cashflows.coupons)
    # the first coupon paying after d is the one we are after
    if ~has_occurred(bond.cashflows[i], d, false)
      # if isa(bond.cashflows[i], Coupon)
      return accrued_amount(bond.cashflows[i], d) / notional(bond, d) * 100.0
    end
  end

  return 0.0
end

type CallableBondArgs{F <: Frequency}
  couponDates::Vector{Date}
  couponAmounts::Vector{Float64}
  redemption::Float64
  redemptionDate::Date
  freq::F
  putCallSchedule::CallabilitySchedule
  callabilityPrices::Vector{Float64}
  callabilityDates::Vector{Date}
end

function CallableBondArgs(bond::CallableFixedRateBond)
  redemption = get_redemption(bond)
  settlement = settlement_date(bond)
  coupDates = Vector{Date}()
  coupAmts = Vector{Float64}()
  for i = 1:length(bond.cashflows.coupons)
    if ~has_occurred(bond.cashflows[i], settlement, false)
      push!(coupDates, date(bond.cashflows[i]))
      push!(coupAmts, amount(bond.cashflows[i]))
    end
  end

  callabilityPrices = Vector{Float64}()
  callabilityDates = Vector{Date}()
  for i in eachindex(bond.putCallSchedule)
    if ~has_occurred(bond.putCallSchedule[i], settlement, false)
      push!(callabilityDates, bond.putCallSchedule[i].date)
      push!(callabilityPrices, bond.putCallSchedule[i].price.amount)

      # calling accrued() forces accrued interest to be zero if future option date
      # is also coupon date, so that dirty price = clean price.  Used here because
      # callability is always applied before coupon in tree engine
      if isa(bond.putCallSchedule[i].price.callType, CleanCall)
        callabilityPrices[end] += accrued(bond, bond.putCallSchedule[i].date)
      end
    end
  end

  return CallableBondArgs(coupDates, coupAmts, amount(redemption), date(redemption), get_frequency(bond), bond.putCallSchedule,
                          callabilityPrices, callabilityDates)
end
