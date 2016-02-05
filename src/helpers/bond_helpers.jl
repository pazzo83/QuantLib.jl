# bond helper functions
using JQuantLib.Instruments, JQuantLib.TermStructures, JQuantLib.PricingEngines

export implied_quote, clean_price, dirty_price, settlement_date

function implied_quote(bond::Bond, yts::YieldTermStructure, pe::PricingEngine, clean_price::Bool)
  settlement_value = calculate(pe, yts, bond)
  return clean_price ? clean_price(bond, settlement_value, settlement_date(bond)) : dirty_price(bond, settlement_value, settlement_date(bond))
end

clean_price(bond::Bond, settlement_value::Float64, settlement_date::Date) = dirty_price(bond, settlement_value, settlement_date) - accrued_amount(bond, settlement_date)

dirty_price(bond::Bond, settlement_value::Float64, settlement_date::Date) = settlement_value * 100.0 / bond.faceAmount # replace with notionals

function settlement_date(bond::Bond, d::Date = Date())
  if d == Date()
    return settings.evaluation_date
  end

  return d + get_settlement_days(bond)
end
