# general swaption engine types and methods
struct NullSwaptionEngine <: PricingEngine end

_calculate!(::NullSwaptionEngine, ::Swaption) = error("Must set valid pricing engine, use update_pricing_engine")

function greater_than_or_equal_to(x::T, y::T) where {T}
  return x >= y
end
