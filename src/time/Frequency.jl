# Frequency
using Dates

abstract type Frequency end

struct NoFrequency <: Frequency end
struct Once <: Frequency end
struct Annual <: Frequency end
struct Semiannual <: Frequency end
struct EveryFourthMonth <: Frequency end
struct Quarterly <: Frequency end
struct Bimonthly <: Frequency end
struct Monthly <: Frequency end
struct EveryFourthWeek <: Frequency end
struct Biweekly <: Frequency end
struct Weekly <: Frequency end
struct Daily <: Frequency end
struct OtherFrequency <: Frequency end

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
