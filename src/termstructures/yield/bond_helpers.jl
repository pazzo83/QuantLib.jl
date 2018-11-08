# Bond Helpers
mutable struct FixedRateBondHelper{DC <: DayCount, P <: PricingEngine} <: BondHelper
  price::Quote
  bond::FixedRateBond{DC, P}
end
value(b::FixedRateBondHelper) = b.price.value
maturity_date(b::FixedRateBondHelper) = maturity_date(b.bond)

# bond helper functions
function implied_quote(bond_h::BondHelper, clean::Bool = true)
  bond = bond_h.bond
  recalculate!(bond)
  settlement_value = bond.settlementValue
  return clean ? clean_price(bond, settlement_value, settlement_date(bond)) : dirty_price(bond, settlement_value, settlement_date(bond))
end

function clone(fixedRateBondHelper::FixedRateBondHelper, ts::TermStructure)
  # first we have to clone the PE
  newPE = clone(fixedRateBondHelper.bond.pricingEngine, ts)
  # then we have to clone the bond
  newBond = clone(fixedRateBondHelper.bond, newPE)

  return FixedRateBondHelper{typeof(newBond.dc), typeof(newPE)}(fixedRateBondHelper.price, newBond)
end

update_termstructure(bondHelper::FixedRateBondHelper, ts::TermStructure) = clone(bondHelper, ts)
