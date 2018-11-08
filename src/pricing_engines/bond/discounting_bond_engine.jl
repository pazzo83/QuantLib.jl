using Dates

struct DiscountingBondEngine{Y <: YieldTermStructure} <: PricingEngine
  yts::Y

  # function call(::Type{DiscountingBondEngine})
  #   new{YieldTermStructure}()
  # end
  #
  # function call{Y}(::Type{DiscountingBondEngine}, yts::Y)
  #   new{Y}(yts)
  # end
end

DiscountingBondEngine() = DiscountingBondEngine{NullYieldTermStructure}(NullYieldTermStructure())

function _calculate!(pe::DiscountingBondEngine, bond::Bond)
  yts = pe.yts
  valuation_date = yts.referenceDate
  value = npv(bond.cashflows, yts, valuation_date, valuation_date)
  bond.settlementValue = value
  # if the valuation_date is the same as the bond settlement date, then we don't have to recalculate

  return bond
end

clone(pe::DiscountingBondEngine, ts::Y) where {Y <: YieldTermStructure} = DiscountingBondEngine{Y}(ts)
