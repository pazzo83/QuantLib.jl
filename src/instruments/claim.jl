using Dates

struct FaceValueClaim <: AbstractClaim end

amount(::FaceValueClaim, ::Date, notional::Float64, recoveryRate::Float64) = notional * (1.0 - recoveryRate)
