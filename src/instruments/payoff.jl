type PlainVanillaPayoff{OT <: OptionType} <: StrikedTypePayoff
  optionType::OT
  strike::Float64
end
