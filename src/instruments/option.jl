type Put <: OptionType end
type Call <: OptionType end

value(::Put) = -1
value(::Call) = 1
