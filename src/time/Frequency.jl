# Frequency
using Base.Dates

abstract Frequency

immutable NoFrequency <: Frequency end
immutable Once <: Frequency end
immutable Annual <: Frequency end
immutable Semiannual <: Frequency end
immutable EveryFourthMonth <: Frequency end
immutable Quarterly <: Frequency end
immutable Bimonthly <: Frequency end
immutable Monthly <: Frequency end
immutable EveryFourthWeek <: Frequency end
immutable Biweekly <: Frequency end
immutable Weekly <: Frequency end
immutable Daily <: Frequency end
immutable OtherFrequency <: Frequency end

# default value if no freq specified
value(::Frequency) = -1

value(::NoFrequency) = -1
value(::Once)			 	= 0
value(::Annual)			= 1
value(::Semiannual)		= 2
value(::EveryFourthMonth)  = 3
value(::Quarterly)		 	= 4
value(::Bimonthly)		 	= 6
value(::Monthly)			= 12
value(::EveryFourthWeek)  	= 13
value(::Biweekly)		 	= 26
value(::Weekly)			= 52
value(::Daily)			 	= 365
value(::OtherFrequency)   	= 999

period(::Annual) = Year(1)
period(::Semiannual) = Month(6)
