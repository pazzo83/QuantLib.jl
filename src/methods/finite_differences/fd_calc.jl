type FdmAffineModelTermStructure{I <: Integer, B <: BusinessCalendar, DC <: DayCount, A <: AffineModel} <: YieldTermStructure
  settlement_days::I
  referenceDate::Date
  calendar::B
  dc::DC
  modelReferenceDate::Date
  model::A
  r::Vector{Float64}
  t::Float64
end

FdmAffineModelTermStructure{B <: BusinessCalendar, DC <: DayCount, A <: AffineModel}(referenceDate::Date, cal::B, dc::DC, modelReferenceDate::Date, model::A, r::Vector{Float64}) =
                            FdmAffineModelTermStructure{Int, B, DC, A}(0, referenceDate, cal, dc, modelReferenceDate, model, r, year_fraction(dc, modelReferenceDate, referenceDate))

discount_impl(ts::FdmAffineModelTermStructure, T::Float64) = discount_bond(ts.model, ts.t, T + ts.t, ts.r)

set_variable(ts::FdmAffineModelTermStructure, r::Vector{Float64}) = ts.r = r

type FdmAffineModelSwapInnerValue{M1 <: Model, M2 <: Model, FM <: FdmMesher, I <: Integer} <: FdmInnerValueCalculator
  disModel::M1
  fwdModel::M2
  swap::VanillaSwap
  exerciseDates::Dict{Float64, Date}
  mesher::FM
  direction::I
  disTs::FdmAffineModelTermStructure
  fwdTs::FdmAffineModelTermStructure

  function FdmAffineModelSwapInnerValue{M1, M2, FM, I}(disModel::M1, fwdModel::M2, swap::VanillaSwap, exerciseDates::Dict{Float64, Date}, mesher::FM, direction::I)
    idx = swap.iborIndex
    newIbor = IborIndex(idx.familyName, idx.tenor, idx.fixingDays, idx.currency, idx.fixingCalendar, idx.convention, idx.endOfMonth, idx.dc)
    newSwap = VanillaSwap(swap.swapT, swap.nominal, swap.fixedSchedule, swap.fixedRate, swap.fixedDayCount, newIbor, swap.spread, swap.floatSchedule,
                          swap.floatDayCount, swap.pricingEngine, swap.paymentConvention)

    return new{M1, M2, FM, I}(disModel, fwdModel, newSwap, exerciseDates, mesher, direction)
  end
end

FdmAffineModelSwapInnerValue{M1 <: Model, M2 <: Model, FM <: FdmMesher, I <: Integer}(disModel::M1, fwdModel::M2, swap::VanillaSwap, exerciseDates::Dict{Float64, Date}, mesher::FM, direction::I) =
                            FdmAffineModelSwapInnerValue{M1, M2, FM, I}(disModel, fwdModel, swap, exerciseDates, mesher, direction)

function get_state{I <: Integer}(::G2, calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, ::Float64)
  retVal = Vector{Float64}(2)
  retVal[1] = get_location(calc.mesher, coords, calc.direction)
  retVal[2] = get_location(calc.mesher, coords, calc.direction + 1)

  return retVal
end

get_state{I <: Integer}(model::HullWhite, calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, t::Float64) = [short_rate(get_dynamics(model), t, get_location(calc.mesher, coords, calc.direction))]

function inner_value{I <: Integer}(calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, i::I, t::Float64)
  iterExerciseDate = get(calc.exerciseDates, t, Date())
  disRate = get_state(calc.disModel, calc, coords, t)
  fwdRate = get_state(calc.fwdModel, calc, coords, t)

  if !isdefined(calc, :disTs) || iterExerciseDate != reference_date(calc.disTs)
    disc = calc.disModel.ts
    calc.disTs = FdmAffineModelTermStructure(iterExerciseDate, disc.calendar, disc.dc, reference_date(disc), calc.disModel, disRate)
    fwd = calc.fwdModel.ts
    calc.fwdTs = FdmAffineModelTermStructure(iterExerciseDate, fwd.calendar, fwd.dc, reference_date(fwd), calc.fwdModel, fwdRate)
    # probably should put this in something more logical
    calc.swap.iborIndex.ts = calc.fwdTs
  else
    # do some updating
    set_variable(calc.disTs, disRate)
    set_variable(calc.fwdTs, fwdRate)
  end

  npv = 0.0
  for j = 1:2
    for i = 1:length(calc.swap.legs[j].coupons)
      cf = calc.swap.legs[j].coupons[i]
      if isa(cf, Coupon)
        npv += accrual_start_date(cf) >= iterExerciseDate ? amount(cf) * discount(calc.disTs, date(cf)) : 0.0
      end
    end
    if j == 1
      npv *= -1.0
    end
  end
  if isa(calc.swap.swapT, Receiver)
    npv *= -1.0
  end
  return max(0.0, npv)
end

avg_inner_value{I <: Integer}(calc::FdmAffineModelSwapInnerValue, coords::Vector{I}, i::I, t::Float64) = inner_value(calc, coords, i, t)
