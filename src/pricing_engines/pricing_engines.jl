# Pricing Engines module
# module PricingEngines
type NullPricingEngine <: PricingEngine end

_calculate!(pe::NullPricingEngine, inst::Instrument) = error("You must set a valid pricing engine")
