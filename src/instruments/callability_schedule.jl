type DirtyCall <: CallType end
type CleanCall <: CallType end

type Price{T <: CallType}
  amount::Float64
  callType::T
end

type Callability{T <: OptionType}
  price::Price
  optionType::T
  date::Date
end

typealias CallabilitySchedule Vector{Callability}
