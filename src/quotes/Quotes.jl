# Quotes
# module Quotes # not a module anymore
type Quote
  value::Float64
  is_valid::Bool

  Quote(value::Float64, is_valid::Bool = true) = new(value, is_valid)
end

# end
