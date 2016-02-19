# Curves

type NullCurve <: Curve end

type FittingCost{I <: Integer} <: CostFunction
  # value::Vector{Float64}
  # values::Vector{Float64}
  firstCashFlow::Vector{I}
  curve::Curve
end

function FittingCost{I <: Integer, C <: Curve}(size::I, curve::C)
  # value = Vector{Float64}()
  # values = Vector{Float64}()
  firstCashFlow = zeros(I, size)

  return FittingCost(firstCashFlow, curve)
end

# Interpolated Curve methods #
max_date{C <: InterpolatedCurve}(curve::C) = curve.dates[end]

function discount{C <: Curve}(curve::C, t::Float64)
  return discount_impl(curve, t)
end

function discount_impl{C <: InterpolatedCurve}(curve::C, t::Float64)
  if t <= curve.times[end]
    return JQuantLib.Math.value(curve.interp, t)
  end

  # println("outside")
  # println(t)
  # println(curve.times)

  # do flat fwd extrapolation
end

function perform_calculations!{C <: InterpolatedCurve}(curve::C)
  _calculate!(curve.boot, curve)
  return curve
end

function value{C <: CostFunction, T}(cf::C, x::Vector{T})
  ref_date = cf.curve.referenceDate
  dc = cf.curve.dc
  squared_error = 0.0
  n = length(cf.curve.bonds)

  # for (i, bh) in enumerate(cf.curve.bonds)
  for i=1:length(cf.curve.bonds)
    bond = cf.curve.bonds[i].bond
    bond_settlement = get_settlement_date(bond)
    model_price = -accrued_amount(bond, bond_settlement)
    leg = bond.cashflows
    for k = cf.firstCashFlow[i]:length(leg.coupons)
      # @inbounds df = discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.coupons[k])))
      @inbounds model_price += amount(leg.coupons[k]) * discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.coupons[k])))
    end

    # redemption
    # @inbounds model_price += amount(leg.redemption) * discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, date(leg.redemption)))

    # adjust NPV for forward settlement
    if bond_settlement != ref_date
      model_price /= discount_function(cf.curve.fittingMethod, x, year_fraction(dc, ref_date, bond_settlement))
    end

    market_price = bond.faceAmount
    price_error = model_price - market_price
    @inbounds weighted_error = cf.curve.fittingMethod.commons.weights[i] * price_error
    squared_error += weighted_error * weighted_error
  end

  return squared_error
end
