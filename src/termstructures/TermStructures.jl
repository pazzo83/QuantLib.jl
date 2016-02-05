# TermStructures module
# module TermStructures
# core TermStructure methods
function check_range{T <: TermStructure}(ts::T, date::Date)
  date < ts.referenceDate && "Date $date before reference_date $(ts.referenceDate)"
end

function check_range{T <: TermStructure}(ts::T, time_frac::Float64)
  time_frac < 0.0 && "Negative time given: $time"
end

max_date{T <: TermStructure}(ts::T) = ts.referenceDate + Date.Year(100)

time_from_reference{T <: TermStructure}(ts::T, date::Date) = year_fraction(ts.dc, reference_date(ts), date, Date(), Date())

function reference_date{T <: TermStructure}(ts::T)
  if ts.referenceDate == Date()
    ts.referenceDate = advance(Dates.Day(ts.settlementDays), ts.calendar, settings.evaluation_date)
  end

  return ts.referenceDate
end

# end
