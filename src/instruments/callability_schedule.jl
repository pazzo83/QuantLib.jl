struct DirtyCall <: CallType end
struct CleanCall <: CallType end

struct Price{T <: CallType}
  amount::Float64
  callType::T
end

struct Callability{T <: OptionType}
  price::Price
  optionType::T
  date::Date
end

const CallabilitySchedule = Vector{Callability}

function has_occurred(callability::Callability, ref_date::Date, include_ref_date::Bool = true)
  # will need to expand this
  if ref_date < callability.date || (ref_date == callability.date && include_ref_date)
    return false
  else
    return true
  end
end
