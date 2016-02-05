# cash_flows.jl
# module CF

using JQuantLib.Time, JQuantLib.Math

# BlackIborCouponPricer() = BlackIborCouponPricer(0.0, 0.0, false)

type CouponMixin{DC <: DayCount}
  accrualStartDate::Date
  accrualEndDate::Date
  refPeriodStart::Date
  refPeriodEnd::Date
  dc::DC
  accrualPeriod::Float64
end

accrual_start_date{C <: Coupon}(coup::C) = coup.couponMixin.accrualStartDate
accrual_end_date{C <: Coupon}(coup::C) = coup.couponMixin.accrualEndDate
ref_period_start{C <: Coupon}(coup::C) = coup.couponMixin.refPeriodStart
ref_period_end{C <: Coupon}(coup::C) = coup.couponMixin.refPeriodEnd
get_dc{C <: Coupon}(coup::C) = coup.couponMixin.dc

# accural_period{C <: Coupon}(coup::C) = coup.couponMixin.accrualPeriod
accrual_period!{C <: Coupon}(coup::C, val::Float64) = coup.couponMixin.accrualPeriod = val
function accrual_period{C <: Coupon}(coup::C)
  if coup.couponMixin.accrualPeriod == -1.0
    p = year_fraction(get_dc(coup), accrual_start_date(coup), accrual_end_date(coup), ref_period_start(coup), ref_period_end(coup))
    accrual_period!(coup, p)
  end

  return coup.couponMixin.accrualPeriod
end

## types of cash flows ##
type SimpleCashFlow <: CashFlow
  amount::Float64
  date::Date
end

amount(cf::SimpleCashFlow) = cf.amount
date(cf::SimpleCashFlow) = cf.date
date_accrual_end(cf::SimpleCashFlow) = cf.date

date{C <: Coupon}(coup::C) = coup.paymentDate
date_accrual_end{C <: Coupon}(coup::C) = accrual_end_date(coup)

type Dividend <: CashFlow
  amount::Float64
  date::Date
end

amount(div::Dividend) = div.amount
date(div::Dividend) = div.date
date_accrual_end(div::Dividend) = div.date

# legs to build cash flows
abstract Leg <: CashFlows

type ZeroCouponLeg <: Leg
  redemption::SimpleCashFlow
end

## Function wrapper for solvers ##
type IRRFinder{L <: Leg, DC <: DayCount, C <: CompoundingType, F <: Frequency}
  leg::L
  npv::Float64
  dc::DC
  comp::C
  freq::F
  includeSettlementDateFlows::Bool
  settlementDate::Date
  npvDate::Date
end

## this function can pass itself or its derivative ##
function operator(finder::IRRFinder)
  function _inner(y::Float64)
    yld = InterestRate(y, finder.dc, finder.comp, finder.freq)
    NPV = npv(finder.leg, yld, finder.includeSettlementDateFlows, finder.settlementDate, finder.npvDate)
    return finder.npv - NPV
  end

  # derivative
  function _inner(y::Float64, ::Derivative)
    yld = InterestRate(y, finder.dc, finder.comp, finder.freq)
    return duration(ModifiedDuration(), finder.leg, yld, finder.dc, finder.includeSettlementDateFlows, finder.settlementDate, finder.npvDate)
  end

  return _inner
end

get_latest_coupon{L <: Leg}(leg::L) = get_latest_coupon(leg, leg.coupons[end])
get_latest_coupon{L <: Leg}(leg::L, simp::SimpleCashFlow) = leg.coupons[end - 1]
get_latest_coupon{L <: Leg, C <: Coupon}(leg::L, coup::C) = coup

check_coupon{C <: CashFlow}(x::C) = isa(x, Coupon)

get_pay_dates{C <: Coupon}(coups::Vector{C}) = Date[date(coup) for coup in coups]

get_reset_dates{C <: Coupon}(coups::Vector{C}) = Date[accrual_start_date(coup) for coup in coups]


## Pricer Methods ##
# function initialize!(pricer::IborCouponPricer, coupon::IborCoupon)
#   # stuff
#   if !pricer.initialized
#     payment_date = date(coupon)
#     if payment_date > coupon.iborIndex.yts.referenceDate
#       pricer.discount = discount(yts, payment_date)
#     else
#       pricer.discount = 1.0
#     end
#     pricer.spreadLegValue = coupon.spread * accrual_period(coupon) * pricer.discount
#
#     pricer.initialized = true
#   end
#
#   return pricer
# end
#
# function adjusted_fixing(pricer::BlackIborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure,
#                         fixing::Float64 = index_fixing(coupon, yts))
#   # stuff
#   if !coupon.isInArrears
#     return fixing
#   end
#
#   d1 = coupon.fixingDate
#   ref_date = cap_vol.referenceDate
#
#   if d1 <= ref_date
#     return fixing
#   end
#
#   d2 = value_date(coupon.iborIndex, d1)
#   d3 = maturity_date(coupon.iborIndex, d2)
#   tau = year_fraction(coupon.iborIndex.dc, d2, d3)
#   varience = black_varience(cap_vol, d1, fixing)
#   adjustment = fixing * fixing * varience * tau / (1.0 + fixing * tau)
#   return fixing + adjustment
# end
#
# function swaplet_price(pricer::BlackIborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure)
#   swapl_price = adjusted_fixing(pricer, coupon, yts, cap_vol) * accrual_period(coupon) * pricer.discount
#   return coupon.gearing * swapl_price + pricer.spreadLegValue
# end
#
# swaplet_rate(pricer::BlackIborCouponPricer, coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure) =
#   swaplet_price(pricer, coupon, yts) / (accrual_period(coupon) * pricer.discount)


## Floating Rate Coupon Methods ##
# function calc_rate(coupon::IborCoupon, yts::YieldTermStructure, cap_vol::OptionletVolatilityStructure)
#   # first initialize pricer
#   initialize!(coupon.pricer, coupon, yts)
#
#   # then get swaplet rate
#   return swaplet_rate(coupon.pricer, coupon, yts, cap_vol)
# end

## NPV METHODS ##
function npv{L <: Leg, Y <: YieldTermStructure}(leg::L, yts::Y, settlement_date::Date, npv_date::Date)
  # stuff
  totalNPV = 0.0
  if length(leg.coupons) == 0
    return totalNPV
  end

  for i in leg
    if has_occurred(i, settlement_date)
      continue
    end
    @inbounds totalNPV += amount(i) * discount(yts, date(i))
  end

  # redemption - not needed with new iterator
  # totalNPV += amount(leg.redemption) * discount(yts, leg.redemption.date)
  return totalNPV / discount(yts, npv_date)
end

function npv{L <: Leg}(leg::L, y::InterestRate, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date)
  if length(leg.coupons) == 0
    return 0.0
  end

  totalNPV = 0.0
  discount_ = 1.0
  last_date = npv_date

  for cp in leg
    if has_occurred(cp, settlement_date)
      continue
    end
    coupon_date = date(cp)
    amount_ = amount(cp)

    if isa(cp, Coupon)
      ref_start_date = ref_period_start(cp)
      ref_end_date = ref_period_end(cp)
    else
      if last_date == npv_date
        ref_start_date = coupon_date - Dates.Year(1)
      else
        ref_start_date = last_date
      end
      ref_end_date = coupon_date
    end

    b = discount_factor(y, last_date, coupon_date, ref_start_date, ref_end_date)

    discount_ *= b
    last_date = coupon_date

    totalNPV += amount_ * discount_
  end

  # redemption - not needed with iterator
  #redempt_date = leg.redemption.date
  # amount = amount(leg.redemption)

  # b = discount_factor(leg.redemption, last_date, redempt_date, last_date, redempt_date)
  # discount *= b
  # totalNPV += amount * discount

  # now we return total npv
  return totalNPV
end

function npv{Y <: YieldTermStructure}(leg::ZeroCouponLeg, yts::Y, settlement_date::Date, npv_date::Date)
  if amount(leg.redemption) == 0.0
    return 0.0
  end

  totalNPV = amount(leg.redemption) * discount(yts, date(leg.redemption))
  return totalNPV / discount(yts, npv_date)
end

function npv(leg::ZeroCouponLeg, y::InterestRate, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date)
  if amount(leg.redemption) == 0.0
    return 0.0
  end

  redempt_date = date(leg.redemption)

  ref_start_date = redempt_date - Dates.Year(1)
  ref_end_date = redempt_date

  discount_ = discount_factor(y, npv_date, redempt_date, ref_start_date, ref_end_date)

  totalNPV = amount(leg.redemption) * discount_

  return totalNPV
end


function npvbps{L <: Leg, Y <: YieldTermStructure}(leg::L, yts::Y, settlement_date::Date, npv_date::Date)
  npv = 0.0
  bps = 0.0

  for cp in leg
    if has_occurred(cp, settlement_date)
      continue
    end

    df = discount(yts, date(cp))
    npv += amount(cp) * df
    if isa(cp, Coupon)
      bps += cp.nominal * accrual_period(cp) * df
    end
    # if isa(cp, FixedRateCoupon)
    #   println("###########################")
    #   println("df ", df)
    #   println("bps ", bps)
    #   println("Date ", date(cp))
    #   println("nom ", cp.nominal)
    #   println("acru p ", accrual_period(cp))
    #   println("DayCount ", cp.dc)
    # end
  end
  d = discount(yts, npv_date)
  npv /= d
  bps = BASIS_POINT * bps / d

  return npv, bps
end

## Duration Calculations ##
modified_duration_calc{F <: Frequency}(::SimpleCompounding, c::Float64, B::Float64, t::Float64, ::Float64, ::F) = c * B * B * t
modified_duration_calc{F <: Frequency}(::CompoundedCompounding, c::Float64, B::Float64, t::Float64, r::Float64, N::F) = c * t * B / (1 + r / JQuantLib.Time.value(N))
modified_duration_calc{F <: Frequency}(::ContinuousCompounding, c::Float64, B::Float64, t::Float64, ::Float64, ::F) = c * B * t
modified_duration_calc{F <: Frequency}(::SimpleThenCompounded, c::Float64, B::Float64, t::Float64, r::Float64, N::F) =
  t <= 1.0 / N ? modified_duration_calc(Simple(), c, B, t, r, N) : modified_duration_calc(CompoundedCompounding(), c, B, t, r, N)

function duration{L <: Leg, DC <: DayCount}(::ModifiedDuration, leg::L, y::InterestRate, dc::DC, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date = Date())
  if length(leg.coupons) == 0 # TODO make this applicable to redemption too
    return 0.0
  end

  if npv_date == Date()
    npv_date = settlement_date
  end

  P = 0.0
  t = 0.0
  dPdy = 0.0
  r = y.rate
  N = y.freq
  last_date = npv_date

  for cp in leg
    if has_occurred(cp, settlement_date)
      continue
    end

    c = amount(cp)
    coupon_date = date(cp)
    if isa(cp, FixedRateCoupon)
      ref_start_date = ref_period_start(cp)
      ref_end_date = ref_period_end(cp)
    else
      if last_date == npv_date
        ref_start_date = coupon_date - Dates.Year(1)
      else
        ref_start_date = last_date
      end
      ref_end_date = coupon_date
    end

    t += year_fraction(dc, last_date, coupon_date, ref_start_date, ref_end_date)

    B = discount_factor(y, t)
    P += c * B

    dPdy -= modified_duration_calc(y.comp, c, B, t, r, N)

    last_date = coupon_date
  end

  if P == 0.0
    return 0.0 # no cashflows
  end

  return -dPdy / P # reverse derivative sign
end

function duration{DC <: DayCount}(::ModifiedDuration, leg::ZeroCouponLeg, y::InterestRate, dc::DC, include_settlement_cf::Bool, settlement_date::Date, npv_date::Date = Date())
  if amount(leg.redemption) == 0.0 || has_occurred(leg.redemption, settlement_date)
    return 0.0
  end

  redempt_date = date(leg.redemption)

  if npv_date == Date()
    npv_date = settlement_date
  end

  r = y.rate
  N = y.freq
  c = amount(leg.redemption)

  ref_start_date = redempt_date - Dates.Year(1)
  ref_end_date = redempt_date

  t = year_fraction(dc, npv_date, redempt_date, ref_start_date, ref_end_date)

  B = discount_factor(y, t)
  P = c * B

  dPdy = modified_duration_calc(y.comp, c, B, t, r, N)

  if P == 0.0
    return 0.0
  end

  return -dPdy / P # reverse derivative sign
end

function aggregate_rate{L <: Leg, I <: Integer}(cf::L, _first::I, _last::I)
  if cf.coupons[_first] == cf.coupons[_last]
    return 0.0
  end

  _next_cf = cf.coupons[_first]
  paymentDate = date(_next_cf)
  result = 0.0
  i = _first
  while i < length(cf.coupons) && date(_next_cf) == paymentDate
    if isa(_next_cf, Coupon)
      result += calc_rate(_next_cf)
    end

    i += 1
    _next_cf = cf.coupons[i]
  end

  return result
end


# functions for sorting and finding
sort_cashflow{C <: Coupon}(cf::C) = accrual_end_date(cf) #cf.couponMixin.accrualEndDate
sort_cashflow(simp::SimpleCashFlow) = simp.date

prev_cf{C <: Coupon}(cf::C, d::Date) = d > cf.paymentDate
prev_cf(simp::SimpleCashFlow) = d > simp.date

next_cf{C <: Coupon}(cf::C, d::Date) = d < cf.paymentDate
next_cf(simp::SimpleCashFlow) = d < simp.date

function previous_cashflow_date{L <: Leg}(cf::L, settlement_date::Date)
  # right now we can assume cashflows are sorted by date because of schedule
  prev_cashflow_idx = findprev(prev_cf, cf.coupons, length(cf.coupons), settlement_date)

  return prev_cashflow_idx == 0 ? 0 : cf.coupons[prev_cashflow_idx].paymentDate
end

function accrual_days{C <: CashFlows, DC <: DayCount}(cf::C, dc::DC, settlement_date::Date)
  last_payment = previous_cashflow_date(cf, settlement_date)

  if last_payment == 0
    return 0.0
  else
    return day_count(dc, last_payment, settlement_date)
  end
end

function next_cashflow{L <: Leg}(cf::L, settlement_date::Date)
  if settlement_date > date(cf.coupons[end])
    return length(cf.coupons)
  end

  return findnext(next_cf, cf.coupons, 1, settlement_date)
end

function next_coupon_rate{L <: Leg}(cf::L, settlement_date::Date)
  next_cf_idx = next_cashflow(cf, settlement_date)
  return aggregate_rate(cf, next_cf_idx, length(cf.coupons))
end

function accrued_amount{L <: Leg}(cf::L, settlement_date::Date, include_settlement_cf::Bool = false)
  next_cf_idx = next_cashflow(cf, settlement_date)
  if cf.coupons[next_cf_idx] == length(cf.coupons)
    return 0.0
  end

  _next_cf = cf.coupons[next_cf_idx]
  paymentDate = date(_next_cf)
  result = 0.0
  i = next_cf_idx
  while i < length(cf.coupons) && date(_next_cf) == paymentDate
    result += accrued_amount(_next_cf, settlement_date)
    i += 1
    _next_cf = cf.coupons[i]
  end

  return result
end

accrued_amount(cf::ZeroCouponLeg, ::Date, ::Bool= false) = 0.0
accrued_amount(simp::SimpleCashFlow, ::Date, ::Bool = false) = 0.0

# accrual_period!{C <: Coupon}(coup::C, p::Float64) = coup.accrualPeriod = p
# function accrual_period{C <: Coupon}(coup::C)
#   if accrual_period(coup) == -1.0
#     p = year_fraction(get_dc(coup), accrual_start_date(coup), accrual_end_date(coup), ref_period_start(coup), ref_period_end(coup))
#     accrual_period!(coup, p)
#   end
#
#   return p
# end

function has_occurred{C <: CashFlow}(cf::C, ref_date::Date, include_settlement_cf::Bool = true)
  # will need to expand this
  if ref_date < date(cf) || (ref_date == date(cf) && include_settlement_cf)
    return false
  else
    return true
  end
end

function maturity_date{L <: Leg}(leg::L)
  # sort cashflows
  sorted_cf = sort(leg.coupons, by=sort_cashflow)
  return date_accrual_end(sorted_cf[end])
end

# Yield Calculations ##
function yield{L <: Leg, DC <: DayCount, C <: CompoundingType, F <: Frequency, I <: Integer}(leg::L, npv::Float64, dc::DC, compounding::C, freq::F, include_settlement_cf::Bool, settlement_date::Date,
              npv_date::Date, accuracy::Float64, max_iter::I, guess::Float64)
  solver = NewtonSolver(max_iter)
  obj_fun = IRRFinder(leg, npv, dc, compounding, freq, include_settlement_cf, settlement_date, npv_date)

  return solve(solver, JQuantLib.operator(obj_fun), accuracy, guess, guess / 10.0)
end

## ITERATORS ##
# this is to iterate through cash flows and redemption
# Base.start(f::FixedRateLeg) = 1
# function Base.next(f::FixedRateLeg, state)
#   if state > length(f.coupons)
#     f.redemption, state + 1
#   else
#     f.coupons[state], state + 1
#   end
# end
#
# Base.done(f::FixedRateLeg, state) = length(f.coupons) + 1 < state

Base.start{L <: Leg}(f::L) = 1
Base.next{L <: Leg}(f::L, state) = f.coupons[state], state + 1
Base.done{L <: Leg}(f::L, state) = length(f.coupons) == state - 1

# end
