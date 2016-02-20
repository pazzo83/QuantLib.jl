using JQuantLib.Time, JQuantLib.Math

const BASIS_POINT = 0.0001

## TYPES ##
## Coupon pricers
type BlackIborCouponPricer{O <: OptionletVolatilityStructure} <: IborCouponPricer
  discount::Float64
  spreadLegValue::Float64
  accrual_period::Float64
  initialized::Bool
  capletVolatility::O

  function call(::Type{BlackIborCouponPricer})
    new{OptionletVolatilityStructure}(0.0, 0.0, 0.0, false)
  end

  function call{O}(::Type{BlackIborCouponPricer}, capletVolatility::O)
    new{O}(0.0, 0.0, 0.0, false, capletVolatility)
  end
end

type IborCoupon{DC <: DayCount, I <: Integer, X <: InterestRateIndex, ICP <: IborCouponPricer} <: Coupon
  couponMixin::CouponMixin{DC}
  paymentDate::Date
  nominal::Float64
  fixingDate::Date
  fixingValueDate::Date
  fixingEndDate::Date
  fixingDays::I
  iborIndex::X
  gearing::Float64
  spread::Float64
  isInArrears::Bool
  spanningTime::Float64
  pricer::ICP
end

function IborCoupon{I <: Integer, X <: InterestRateIndex, DC <: DayCount, ICP <: IborCouponPricer}(paymentDate::Date, nominal::Float64, startDate::Date, endDate::Date, fixingDays::I, iborIndex::X,
                    gearing::Float64, spread::Float64, refPeriodStart::Date, refPeriodEnd::Date, dc::DC, isInArrears::Bool, pricer::ICP)
  # TODO check if right fixing days
  _fixing_date = isInArrears ? fixing_date(iborIndex, endDate) : fixing_date(iborIndex, startDate)
  fixing_cal = iborIndex.fixingCalendar
  idx_fixing_days = iborIndex.fixingDays
  fixing_val_date = advance(Base.Dates.Day(idx_fixing_days), fixing_cal, _fixing_date, iborIndex.convention)

  if isInArrears
    fixing_end_date = maturity_date(iborIndex, fixing_val_date)
  else
    next_fixing = advance(-Base.Dates.Day(fixingDays), fixing_cal, endDate, iborIndex.convention)
    fixing_end_date = advance(Base.Dates.Day(idx_fixing_days), fixing_cal, next_fixing, iborIndex.convention)
  end

  spanning_time = year_fraction(iborIndex.dc, fixing_val_date, fixing_end_date)

  ## TODO ensure positive (> 0) spanning_time

  return IborCoupon{DC, I, X, ICP}(CouponMixin{DC}(startDate, endDate, refPeriodStart, refPeriodEnd, dc, -1.0), paymentDate, nominal, _fixing_date, fixing_val_date, fixing_end_date, fixingDays, iborIndex, gearing, spread,
                    isInArrears, spanning_time, pricer)
end

# Coupon methods
amount(coup::IborCoupon) = calc_rate(coup) * accrual_period(coup) * coup.nominal

get_pay_dates(coups::Vector{Union{IborCoupon, SimpleCashFlow}}) = Date[date(coup) for coup in filter(check_coupon, coups)]
get_reset_dates(coups::Vector{Union{IborCoupon, SimpleCashFlow}}) = Date[accrual_start_date(coup) for coup in filter(check_coupon, coups)]
get_gearings(coups::Vector{Union{IborCoupon, SimpleCashFlow}}) = Float64[coup.gearing for coup in filter(check_coupon, coups)]

function calc_rate(coup::IborCoupon)
  JQuantLib.initialize!(coup.pricer, coup)
  # if coup.isInArrears
  #   println("sr ", sr)
  #   error("BREAK")
  # end
  return swaplet_rate(coup.pricer, coup)
end

function index_fixing(coupon::IborCoupon)
  today = settings.evaluation_date

  if coupon.fixingDate > today
    return forecast_fixing(coupon.iborIndex, coupon.iborIndex.ts, coupon.fixingValueDate, coupon.fixingEndDate, coupon.spanningTime)
  end

  error("Fixing date on or before eval date")
end

function accrued_amount(coup::IborCoupon, settlement_date::Date)
  if settlement_date <= accrual_start_date(coup) || settlement_date > coup.paymentDate
    return 0.0
  end

  return coup.nominal * calc_rate(coup) * year_fraction(get_dc(coup), accrual_start_date(coup), min(settlement_date, accrual_end_date(coup)), ref_period_start(coup), ref_period_end(coup))
end

# Floating legs
type IborLeg <: Leg
  coupons::Vector{Union{IborCoupon, SimpleCashFlow}}

  function IborLeg{X <: InterestRateIndex, DC <: DayCount, C <: BusinessDayConvention, I <: Integer, ICP <: IborCouponPricer}(schedule::Schedule, nominal::Float64, iborIndex::X, paymentDC::DC, paymentAdj::C,
                   fixingDays::Vector{I} = fill(iborIndex.fixingDays, length(schedule.dates) - 1),
                   gearings::Vector{Float64} = ones(length(schedule.dates) - 1), spreads::Vector{Float64} = zeros(length(schedule.dates) - 1),
                   caps::Vector{Float64} = Vector{Float64}(), floors::Vector{Float64} = Vector{Float64}(), isInArrears::Bool = false,
                   isZero::Bool = false, pricer::ICP = BlackIborCouponPricer(); add_redemption::Bool = true)
    n = add_redemption ? length(schedule.dates) : length(schedule.dates) - 1
    coups = Vector{Union{IborCoupon, SimpleCashFlow}}(n)
    last_payment_date = adjust(schedule.cal, paymentAdj, schedule.dates[end])

    _start = ref_start = schedule.dates[1]
    _end = ref_end = schedule.dates[2]
    payment_date = adjust(schedule.cal, paymentAdj, _end)
    #TODO: setup payment adjustments and the like
    coups[1] = IborCoupon(payment_date, nominal, _start, _end, fixingDays[1], iborIndex, gearings[1], spreads[1], ref_start, ref_end,
                          paymentDC, isInArrears, pricer)

    # ref_date = _start = end_date = _end = Date()
    # build coupons
    count = 2
    ref_start = _start = _end
    ref_end = _end = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
    payment_date = adjust(schedule.cal, paymentAdj, _end)

    while _start < schedule.dates[end]
      @inbounds coups[count] = IborCoupon(payment_date, nominal, _start, _end, fixingDays[count], iborIndex, gearings[count], spreads[count], ref_start, ref_end,
                            paymentDC, isInArrears, pricer)

      count += 1
      ref_start = _start = _end
      ref_end = _end = count == length(schedule.dates) ? schedule.dates[end] : schedule.dates[count + 1]
      payment_date = adjust(schedule.cal, paymentAdj, _end)
    end

    # for i = 1:n
    #   ref_start = _start = schedule.dates[i]
    #   ref_end = _end = schedule.dates[i + 1]
    #   payment_date = isZero ? last_payment_date : adjust(schedule.cal, paymentAdj, _end)
    #
    #   # Just Floating Rate Coupons right now
    #   coups[i] = IborCoupon(payment_date, nominal, _start, _end, fixingDays[i], iborIndex, gearings[i], spreads[i], ref_start, ref_end,
    #                         paymentDC, isInArrears, pricer)
    # end

    if add_redemption
      @inbounds coups[end] = SimpleCashFlow(nominal, _end)
    end

    new(coups)
  end
end

## Pricer Methods ##
function initialize!(pricer::BlackIborCouponPricer, coup::IborCoupon)
  idx = coup.iborIndex
  yts = idx.ts

  payment_date = date(coup)
  if payment_date > yts.referenceDate
    pricer.discount = discount(yts, payment_date)
  else
    pricer.discount = 1.0
  end

  pricer.accrual_period = accrual_period(coup)

  pricer.spreadLegValue = coup.spread * pricer.accrual_period * pricer.discount

  return pricer
end

function swaplet_price(pricer::BlackIborCouponPricer, coup::IborCoupon)
  _swaplet_price = adjusted_fixing(pricer, coup) * pricer.accrual_period * pricer.discount
  # if coup.isInArrears
  #   println("gearing ", coup.gearing)
  #   error("BREAK")
  # end
  return coup.gearing * _swaplet_price + pricer.spreadLegValue
end

swaplet_rate(pricer::BlackIborCouponPricer, coup::IborCoupon) = swaplet_price(pricer, coup) / (pricer.accrual_period * pricer.discount)

function adjusted_fixing(pricer::BlackIborCouponPricer, coup::IborCoupon, fixing::Float64 = -1.0)
  if fixing == -1.0
    fixing = index_fixing(coup)
  end

  if !coup.isInArrears
    return fixing
  end

  d1 = coup.fixingDate
  ref_date = pricer.capletVolatility.referenceDate
  if d1 <= ref_date
    return fixing
  end

  # See Hull, 4th ed., page 550
  idx = coup.iborIndex
  d2 = value_date(idx, d1)
  d3 = maturity_date(idx, d2)
  tau = year_fraction(idx.dc, d2, d3)
  variance = black_variance(pricer.capletVolatility, d1, fixing)
  # println("fixing ", fixing)
  # println("d1 ", d1)
  # println("d2 ", d2)
  # println("d3 ", d3)
  # println("tau ", tau)
  # println("variance ", variance)

  adj = fixing * fixing * variance * tau / (1.0 + fixing * tau)
  # println("adj ", adj)
  # error("BREAK")
  return fixing + adj
end

function update_pricer!{O <: OptionletVolatilityStructure}(leg::IborLeg, opt::O)
  for coup in leg.coupons
    if isa(coup, IborCoupon)
      coup.pricer.capletVolatility = opt
    end
  end

  return leg
end
