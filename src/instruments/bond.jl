type DividendSchedule
  dividends::Vector{Dividend}
end

DividendSchedule() = DividendSchedule(Vector{Dividend}(0))

type BondMixin{I <: Integer}
  settlementDays::I
  # schedule::Schedule
  # cashflows::L
  # dc::DC
  issueDate::Date
  maturityDate::Date
  # pricingEngine::P
end

get_settlement_days{B <: Bond}(bond::B) = bond.bondMixin.settlementDays
get_issue_date{B <: Bond}(bond::B) = bond.bondMixin.issueDate
get_maturity_date{B <: Bond}(bond::B) = bond.bondMixin.maturityDate

type FixedRateBond{I <: Integer, DC <: DayCount, P <: PricingEngine} <: Bond
  lazyMixin::LazyMixin
  bondMixin::BondMixin{I}
  faceAmount::Float64
  schedule::Schedule
  cashflows::FixedRateLeg
  dc::DC
  redemption::Float64
  startDate::Date
  pricingEngine::P
  settlementValue::Float64
end

function FixedRateBond{I <: Integer, DC <: DayCount, B <: BusinessDayConvention, C <: BusinessCalendar, P <: PricingEngine}(settlementDays::I, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DC, paymentConvention::B,
  redemption::Float64, issueDate::Date, calendar::C, pricing_engine::P)
  maturityDate = schedule.dates[end]

  # num_payments = length(schedule.dates)
  # coups = zeros(num_payments + 1)
  #
  # # build coupons
  # for i = 1:num_payments
  #   fixed_leg = FixedRateCoupon(schedule, faceAmount, coup_rate, calendar, paymentConvention)
  #   coups[i] = fixed_leg
  # end

  # add redemption
  # coups[end] = SimpleCashFlow(redemption, maturityDate)

  coups = FixedRateLeg(schedule, faceAmount, coup_rate, calendar, paymentConvention, dc)

  return FixedRateBond{I, DC, P}(LazyMixin(), BondMixin{I}(settlementDays, issueDate, maturityDate), faceAmount, schedule, coups, dc, redemption, schedule.dates[1], pricing_engine, 0.0)
end

type FloatingRateBond{I <: Integer, X <: InterestRateIndex, DC <: DayCount, P <: PricingEngine} <: Bond
  lazyMixin::LazyMixin
  bondMixin::BondMixin{I}
  faceAmount::Float64
  schedule::Schedule
  cashflows::IborLeg
  iborIndex::X
  dc::DC
  fixingDays::I
  gearings::Vector{Float64}
  spreads::Vector{Float64}
  caps::Vector{Float64}
  floors::Vector{Float64}
  inArrears::Bool
  redemption::Float64
  pricingEngine::P
  settlementValue::Float64
end

function FloatingRateBond{I <: Integer, X <: InterestRateIndex, DC <: DayCount, B <: BusinessDayConvention, P <: PricingEngine}(settlementDays::I, faceAmount::Float64, schedule::Schedule, iborIndex::X, dc::DC,
                          convention::B, fixingDays::I, issueDate::Date, pricingEngine::P, inArrears::Bool = false, redemption::Float64 = 100.0,
                          gearings::Vector{Float64} = ones(length(schedule.dates) - 1), spreads::Vector{Float64} = zeros(length(schedule.dates) - 1),
                          caps::Vector{Float64} = Vector{Float64}(), floors::Vector{Float64} = Vector{Float64}())
  maturityDate = schedule.dates[end]
  fixingDaysVect = fill(fixingDays, length(schedule.dates) - 1)

  coups = IborLeg(schedule, faceAmount, iborIndex, dc, convention, fixingDaysVect, gearings, spreads, caps, floors, inArrears; add_redemption=true)
  return FloatingRateBond{I, X, DC, P}(LazyMixin(), BondMixin{I}(settlementDays, issueDate, maturityDate), faceAmount, schedule, coups, iborIndex, dc, fixingDays, gearings, spreads, caps, floors,
                          inArrears, redemption, pricingEngine, 0.0)
end

type ZeroCouponBond{I <: Integer, P <: PricingEngine} <: Bond
  lazyMixin::LazyMixin
  bondMixin::BondMixin{I}
  faceAmount::Float64
  redemption::Float64
  cashflows::ZeroCouponLeg
  settlementValue::Float64
  pricingEngine::P
end

function ZeroCouponBond{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, P <: PricingEngine}(settlementDays::I, calendar::B, faceAmount::Float64, maturityDate::Date, paymentConvention::C=Following(),
                        redemption::Float64=100.0, issueDate::Date=Date(), pe::P = DiscountingBondEngine())
  # build redemption CashFlow
  redemption_cf = ZeroCouponLeg(SimpleCashFlow(redemption, maturityDate))
  return ZeroCouponBond{I, P}(LazyMixin(), BondMixin{I}(settlementDays, issueDate, maturityDate), faceAmount, redemption, redemption_cf, 0.0, pe)
end

get_settlement_date{B <: Bond}(b::B) = get_issue_date(b)

function notional{B <: Bond}(bond::B, d::Date)
  if d > get_maturity_date(bond)
    return 0.0
  else
    return bond.faceAmount
  end
end

# Calculation method
function perform_calculations!{B <: Bond}(bond::B)
  bond.settlementValue = 0.0 # reset - TODO this will be expanded
  _calculate!(bond.pricingEngine, bond)

  return bond
end

# bond methods
accrued_amount{B <: Bond}(bond::B, settlement::Date) = accrued_amount(bond.cashflows, settlement, false) * 100.0 / notional(bond, settlement)

maturity_date{B <: Bond}(bond::B) = maturity_date(bond.cashflows)

function yield{B <: Bond, DC <: DayCount, C <: CompoundingType, F <: Frequency, I <: Integer}(bond::B, clean_price::Float64, dc::DC, compounding::C, freq::F, settlement::Date, accuracy::Float64 = 1.0e-10,
              max_iter::I = 100, guess::Float64 = 0.05)
  dirty_price = clean_price + accrued_amount(bond, settlement)
  dirty_price /= 100.0 / notional(bond, settlement)

  return yield(bond.cashflows, dirty_price, dc, compounding, freq, false, settlement, settlement, accuracy, max_iter, guess)
end

function duration{B <: Bond, D <: Duration, DC <: DayCount}(bond::B, yld::InterestRate, duration_::D, dc::DC, settlement_date::Date)
  return duration(duration_, bond.cashflows, yld, dc, false, settlement_date)
end

function duration{B <: Bond, DC <: DayCount, C <: CompoundingType, F <: Frequency, D <: Duration}(bond::B, yld::Float64, dc::DC, compounding::C, freq::F, duration_::D, settlement_date::Date)
  y = InterestRate(yld, dc, compounding, freq)
  return duration(bond, y, duration_, dc, settlement_date)
end

# function npv{B <: Bond, P <: PricingEngine}(bond::B, pe::P)
#   calculate!(pe, bond)
#   return bond.settlementValue
# end
function npv{B <: Bond}(bond::B)
  calculate!(bond)

  return bond.settlementValue
end

clean_price{B <: Bond}(bond::B, settlement_value::Float64, settlement_date::Date) = dirty_price(bond, settlement_value, settlement_date) - accrued_amount(bond, settlement_date)
dirty_price{B <: Bond}(bond::B, settlement_value::Float64, settlement_date::Date) = settlement_value * 100.0 / bond.faceAmount # replace with notionals

clean_price{B <: Bond}(bond::B) = clean_price(bond, bond.settlementValue, settlement_date(bond))
dirty_price{B <: Bond}(bond::B) = dirty_price(bond, bond.settlementValue, settlement_date(bond))

function settlement_date{B <: Bond}(bond::B, d::Date = Date())
  if d == Date()
    d = settings.evaluation_date
  end

  return d + Base.Dates.Day(get_settlement_days(bond))
end
