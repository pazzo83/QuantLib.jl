type DividendSchedule
  dividends::Vector{Dividend}
end

DividendSchedule() = DividendSchedule(Vector{Dividend}(0))

type BondMixin
  settlementDays::Int
  # schedule::Schedule
  # cashflows::L
  # dc::DC
  issueDate::Date
  maturityDate::Date
  # pricingEngine::P
end

get_settlement_days(bond::Bond) = bond.bondMixin.settlementDays
get_issue_date(bond::Bond) = bond.bondMixin.issueDate
get_maturity_date(bond::Bond) = bond.bondMixin.maturityDate

type FixedRateBond{DC <: DayCount, P <: PricingEngine} <: Bond
  lazyMixin::LazyMixin
  bondMixin::BondMixin
  faceAmount::Float64
  schedule::Schedule
  cashflows::FixedRateLeg
  dc::DC
  redemption::Float64
  startDate::Date
  pricingEngine::P
  settlementValue::Float64
end

function FixedRateBond{DC <: DayCount, B <: BusinessDayConvention, C <: BusinessCalendar, P <: PricingEngine}(settlementDays::Int, faceAmount::Float64, schedule::Schedule, coup_rate::Float64, dc::DC, paymentConvention::B,
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

  return FixedRateBond{DC, P}(LazyMixin(), BondMixin(settlementDays, issueDate, maturityDate), faceAmount, schedule, coups, dc, redemption, schedule.dates[1], pricing_engine, 0.0)
end

type FloatingRateBond{X <: InterestRateIndex, DC <: DayCount, P <: PricingEngine} <: Bond
  lazyMixin::LazyMixin
  bondMixin::BondMixin
  faceAmount::Float64
  schedule::Schedule
  cashflows::IborLeg
  iborIndex::X
  dc::DC
  fixingDays::Int
  gearings::Vector{Float64}
  spreads::Vector{Float64}
  caps::Vector{Float64}
  floors::Vector{Float64}
  inArrears::Bool
  redemption::Float64
  pricingEngine::P
  settlementValue::Float64
end

function FloatingRateBond{X <: InterestRateIndex, DC <: DayCount, B <: BusinessDayConvention, P <: PricingEngine}(settlementDays::Int, faceAmount::Float64, schedule::Schedule, iborIndex::X, dc::DC,
                          convention::B, fixingDays::Int, issueDate::Date, pricingEngine::P, inArrears::Bool = false, redemption::Float64 = 100.0,
                          gearings::Vector{Float64} = ones(length(schedule.dates) - 1), spreads::Vector{Float64} = zeros(length(schedule.dates) - 1),
                          caps::Vector{Float64} = Vector{Float64}(), floors::Vector{Float64} = Vector{Float64}();
                          cap_vol::OptionletVolatilityStructure=NullOptionletVolatilityStructure())
  maturityDate = schedule.dates[end]
  fixingDaysVect = fill(fixingDays, length(schedule.dates) - 1)

  coups = IborLeg(schedule, faceAmount, iborIndex, dc, convention, fixingDaysVect, gearings, spreads, caps, floors, inArrears;
                  add_redemption=true, cap_vol = cap_vol)
  return FloatingRateBond{X, DC, P}(LazyMixin(), BondMixin(settlementDays, issueDate, maturityDate), faceAmount, schedule, coups, iborIndex, dc, fixingDays, gearings, spreads, caps, floors,
                          inArrears, redemption, pricingEngine, 0.0)
end

type ZeroCouponBond{BC <: BusinessCalendar, P <: PricingEngine} <: Bond
  lazyMixin::LazyMixin
  bondMixin::BondMixin
  faceAmount::Float64
  redemption::Float64
  cashflows::ZeroCouponLeg
  calendar::BC
  settlementValue::Float64
  pricingEngine::P
end

function ZeroCouponBond{B <: BusinessCalendar, C <: BusinessDayConvention, P <: PricingEngine}(settlementDays::Int, calendar::B, faceAmount::Float64, maturityDate::Date,
                        paymentConvention::C=Following(), redemption::Float64=100.0, issueDate::Date=Date(), pe::P = DiscountingBondEngine())
  # build redemption CashFlow
  redemption_cf = ZeroCouponLeg(SimpleCashFlow(redemption, maturityDate))
  return ZeroCouponBond{B, P}(LazyMixin(), BondMixin(settlementDays, issueDate, maturityDate), faceAmount, redemption, redemption_cf, calendar, 0.0, pe)
end

get_calendar(zcb::ZeroCouponBond) = zcb.calendar

get_calendar(b::Bond) = b.schedule.cal

get_settlement_date(b::Bond) = get_issue_date(b)

function notional(bond::Bond, d::Date)
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
accrued_amount(bond::Bond, settlement::Date) = accrued_amount(bond.cashflows, settlement, false) * 100.0 / notional(bond, settlement)

maturity_date(bond::Bond) = maturity_date(bond.cashflows)

function yield(bond::Bond, clean_price::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, settlement::Date, accuracy::Float64 = 1.0e-10,
              max_iter::Int = 100, guess::Float64 = 0.05)
  dirty_price = clean_price + accrued_amount(bond, settlement)
  dirty_price /= 100.0 / notional(bond, settlement)

  return yield(bond.cashflows, dirty_price, dc, compounding, freq, false, settlement, settlement, accuracy, max_iter, guess)
end

yield(bond::Bond, dc::DayCount, compounding::CompoundingType, freq::Frequency, settlement::Date, accuracy::Float64 = 1.0e-10, max_iter::Int = 100) =
      yield(bond, clean_price(bond), dc, compounding, freq, settlement, accuracy, max_iter)

yield(bond::Bond, dc::DayCount, compounding::CompoundingType, freq::Frequency, accuracy::Float64 = 1.0e-10, max_iter::Int = 100) =
      yield(bond, clean_price(bond), dc, compounding, freq, settlement_date(bond), accuracy, max_iter)

function duration(bond::Bond, yld::InterestRate, duration_::Duration, dc::DayCount, settlement_date::Date)
  return duration(duration_, bond.cashflows, yld, dc, false, settlement_date)
end

function duration(bond::Bond, yld::Float64, dc::DayCount, compounding::CompoundingType, freq::Frequency, duration_::Duration, settlement_date::Date)
  y = InterestRate(yld, dc, compounding, freq)
  return duration(bond, y, duration_, dc, settlement_date)
end

# function npv{B <: Bond, P <: PricingEngine}(bond::B, pe::P)
#   calculate!(pe, bond)
#   return bond.settlementValue
# end
function npv(bond::Bond)
  calculate!(bond)

  return bond.settlementValue
end

clean_price(bond::Bond, settlement_value::Float64, settlement_date::Date) = dirty_price(bond, settlement_value, settlement_date) - accrued_amount(bond, settlement_date)
dirty_price(bond::Bond, settlement_value::Float64, settlement_date::Date) = settlement_value * 100.0 / bond.faceAmount # replace with notionals

function clean_price(bond::Bond)
  calculate!(bond)
  clean_price(bond, bond.settlementValue, settlement_date(bond))
end

dirty_price(bond::Bond) = dirty_price(bond, bond.settlementValue, settlement_date(bond))

function settlement_date(bond::Bond, d::Date = Date())
  if d == Date()
    d = settings.evaluation_date
  end

  # return d + Base.Dates.Day(get_settlement_days(bond))
  return advance(Dates.Day(get_settlement_days(bond)), get_calendar(bond), d)
end

get_redemption(b::Bond) = b.cashflows.coupons[end]
get_frequency(b::Bond) = b.schedule.tenor.freq

## clone methods ##
function clone(bond::FixedRateBond, pe::PricingEngine = bond.pricingEngine)
  # if everything is the same, we keep the calculation status etc
  lazyMixin, settlementValue = pe == bond.pricingEngine ? (bond.lazyMixin, bond.settlementValue) : (LazyMixin(), 0.0)

  return FixedRateBond(lazyMixin, bond.bondMixin, bond.faceAmount, bond.schedule, bond.cashflows, bond.dc, bond.redemption, bond.startDate, pe, settlementValue)
end
