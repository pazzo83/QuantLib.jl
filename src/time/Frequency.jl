# Frequency
using Base.Dates

abstract Frequency

type NoFrequency <: Frequency end
type Once <: Frequency end
type Annual <: Frequency end
type Semiannual <: Frequency end
type EveryFourthMonth <: Frequency end
type Quarterly <: Frequency end
type Bimonthly <: Frequency end
type Monthly <: Frequency end
type EveryFourthWeek <: Frequency end
type Biweekly <: Frequency end
type Weekly <: Frequency end
type Daily <: Frequency end
type OtherFrequency <: Frequency end

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
