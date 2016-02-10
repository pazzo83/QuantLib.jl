type Payer <: SwapType end
type Receiver <: SwapType end

type Buyer <: CDSProtectionSide end

type SwapResults <: Results
  legNPV::Vector{Float64}
  legBPS::Vector{Float64}
  npvDateDiscount::Float64
  startDiscounts::Vector{Float64}
  endDiscounts::Vector{Float64}
  fairRate::Float64
  value::Float64

  function SwapResults{I <: Integer}(n::I)
    legNPV = zeros(n)
    legBPS = zeros(n)
    startDiscounts = zeros(n)
    endDiscounts = zeros(n)

    new(legNPV, legBPS, 0.0, startDiscounts, endDiscounts, -1.0, 0.0)
  end
end

type VanillaSwapArgs
  fixedResetDates::Vector{Date}
  fixedPayDates::Vector{Date}
  floatingResetDates::Vector{Date}
  floatingPayDates::Vector{Date}
  floatingAccrualTimes::Vector{Float64}
  floatingSpreads::Vector{Float64}
  fixedCoupons::Vector{Float64}
  floatingCoupons::Vector{Float64}
end

function VanillaSwapArgs{L <: Leg}(legs::Vector{L})
  fixedCoups = legs[1].coupons
  floatingCoups = legs[2].coupons
  fixedCoupons = [amount(coup) for coup in fixedCoups]
  # floatingCoupons = [amount(coup) for coup in floatingCoups]
  floatingCoupons = zeros(length(floatingCoups))
  floatingAccrualTimes = [accrual_period(coup) for coup in floatingCoups]
  floatingSpreads = [coup.spread for coup in floatingCoups]
  return VanillaSwapArgs(get_reset_dates(fixedCoups), get_pay_dates(fixedCoups), get_reset_dates(floatingCoups), get_pay_dates(floatingCoups), floatingAccrualTimes, floatingSpreads, fixedCoupons, floatingCoupons)
end

function reset!(sr::SwapResults)
  n = length(sr.legNPV)
  sr.legNPV = zeros(n)
  sr.legBPS = zeros(n)
  sr.npvDateDiscount = 0.0
  sr.startDiscounts = zeros(n)
  sr.endDiscounts = zeros(n)
  sr.value = 0.0
  sr.fairRate = 0.0

  return sr
end

type VanillaSwap{ST <: SwapType, DC_fix <: DayCount, DC_float <: DayCount, B <: BusinessDayConvention, L <: Leg, P <: PricingEngine} <: Swap
  lazyMixin::LazyMixin
  swapT::ST
  nominal::Float64
  fixedSchedule::Schedule
  fixedRate::Float64
  fixedDayCount::DC_fix
  iborIndex::IborIndex
  spread::Float64
  floatSchedule::Schedule
  floatDayCount::DC_float
  paymentConvention::B
  legs::Vector{L}
  payer::Vector{Float64}
  pricingEngine::P
  results::SwapResults
  args::VanillaSwapArgs
end

# Constructors
function VanillaSwap{ST <: SwapType, DC_fix <: DayCount, DC_float <: DayCount, B <: BusinessDayConvention, P <: PricingEngine}(swapT::ST, nominal::Float64, fixedSchedule::Schedule, fixedRate::Float64,
                    fixedDayCount::DC_fix, iborIndex::IborIndex, spread::Float64, floatSchedule::Schedule, floatDayCount::DC_float, pricingEngine::P, paymentConvention::B = floatSchedule.convention)
  # build swap cashflows
  legs = Vector{Leg}(2)
  # first leg is fixed
  legs[1] = FixedRateLeg(fixedSchedule, nominal, fixedRate, fixedSchedule.cal, paymentConvention, fixedDayCount; add_redemption=false)
  # second leg is floating
  legs[2] = IborLeg(floatSchedule, nominal, iborIndex, floatDayCount, paymentConvention; add_redemption=false)

  payer = _build_payer(swapT)

  results = SwapResults(2)

  return VanillaSwap{ST, DC_fix, DC_float, B, Leg, P}(LazyMixin(), swapT, nominal, fixedSchedule, fixedRate, fixedDayCount, iborIndex, spread, floatSchedule, floatDayCount,
                    paymentConvention, legs, payer, pricingEngine, results, VanillaSwapArgs(legs))
end

# CDS #
type CreditDefaultSwap{S <: CDSProtectionSide, B <: BusinessDayConvention, DC <: DayCount, P <: PricingEngine} <: Swap
  lazyMixin::LazyMixin
  side::S
  notional::Float64
  spread::Float64
  schedule::Schedule
  convention::B
  dc::DC
  settlesAccrual::Bool
  paysAtDefaultTime::Bool
  protectionStart::Date
  pricingEngine::P

  CreditDefaultSwap(lazyMixin::LazyMixin,
                    side::S,
                    notional::Float64,
                    spread::Float64,
                    schedule::Schedule,
                    convention::B,
                    dc::DC,
                    settlesAccrual::Bool,
                    paysAtDefaultTime::Bool,
                    protectionStart::Date,
                    pricingEngine::P) = new(lazyMixin, side, notional, spread, schedule, convention, dc, settlesAccrual, paysAtDefaultTime, protectionStart, pricingEngine)
end

function CreditDefaultSwap{S <: CDSProtectionSide, B <: BusinessDayConvention, DC <: DayCount, P <: PricingEngine}(side::S, notional::Float64, spread::Float64, schedule::Schedule,
                          convention::B, dc::DC, settlesAccrual::Bool, paysAtDefaultTime::Bool, protectionStart::Date, pricingEngine::P)
  return CreditDefaultSwap{S, B, DC, P}(LazyMixin(), side, notional, spread, schedule, convention, dc, settlesAccrual, paysAtDefaultTime, protectionStart, pricingEngine)
end

# Swap Helper methods
function _build_payer(swapT::Payer)
  x = ones(2)
  x[1] = -1.0
  return x
end

function _build_payer(swapT::Receiver)
  x = ones(2)
  x[2] = -1.0
  return x
end

# Swap methods #
function maturity_date{S <: Swap}(swap::S)
  d = maturity_date(swap.legs[1])
  for i = 2:length(swap.legs)
    d = max(d, maturity_date(swap.legs[i]))
  end

  return d
end

# Calculation method #
function perform_calculations!(swap::VanillaSwap)
  reset!(swap.results) # reset - TODO this will be expanded
  swap.args.floatingCoupons = [amount(coup) for coup in swap.legs[2].coupons]
  _calculate!(swap.pricingEngine, swap)

  return swap
end

floating_leg_NPV(swap::VanillaSwap) = swap.results.legNPV[2]
floating_leg_BPS(swap::VanillaSwap) = swap.results.legBPS[2]

fixed_leg_BPS(swap::VanillaSwap) = swap.results.legBPS[1]

function fair_rate(swap::VanillaSwap)
  calculate!(swap)

  return swap.results.fairRate
end

# some helper methods #
# function clone(swap::VanillaSwap, pe::PricingEngine = swap.pricingEngine)
#   lazyMixin, res, args = pe == swap.pricingEngine ? (swap.lazyMixin, swap.results, swap.args) : (LazyMixin(), SwapResults(2), VanillaSwapArgs(swap.legs))
#   # args = VanillaSwapArgs(swap.legs)
#   # res = SwapResults(2)
#
#   return VanillaSwap(lazyMixin, swap.swapT, swap.nominal, swap.fixedSchedule, swap.fixedRate, swap.fixedDayCount, swap.iborIndex, swap.spread,
#                     swap.floatSchedule, swap.floatDayCount, swap.paymentConvention, swap.legs, swap.payer, pe, res, args)
# end

function clone(swap::VanillaSwap, pe::PricingEngine = swap.pricingEngine, ts::TermStructure = swap.iborIndex.ts)
  is_new = pe != swap.pricingEngine || ts != swap.iborIndex.ts

  lazyMixin, res, args = is_new ? (swap.lazyMixin, swap.results, swap.args) : (LazyMixin(), SwapResults(2), VanillaSwapArgs(swap.legs))

  # we need a new ibor and to rebuild floating rate coupons
  if ts != swap.iborIndex.ts
    newIbor = clone(swap.iborIndex, ts)
    newLegs = Vector{Leg}(2)
    newLegs[1] = swap.legs[1]
    newLegs[2] = IborLeg(swap.floatSchedule, swap.nominal, newIbor, swap.floatDayCount, swap.paymentConvention; add_redemption=false)
  else
    newIbor = swap.iborIndex
    newLegs = swap.legs
  end

  return VanillaSwap(lazyMixin, swap.swapT, swap.nominal, swap.fixedSchedule, swap.fixedRate, swap.fixedDayCount, newIbor, swap.spread,
                    swap.floatSchedule, swap.floatDayCount, swap.paymentConvention, newLegs, swap.payer, pe, res, args)
end

get_pricing_engine_type{ST, DC_fix, DC_float, B, L, P}(::VanillaSwap{ST, DC_fix, DC_float, B, L, P}) = P
