immutable IborIndex{TP <: TenorPeriod, CUR <: AbstractCurrency, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure} <: InterestRateIndex
  familyName::String
  tenor::TP
  fixingDays::Int
  currency::CUR
  fixingCalendar::B
  convention::C
  endOfMonth::Bool
  dc::DC
  ts::T
  pastFixings::Dict{Date, Float64}

  # call{S, I, B, C, DC}(::Type{IborIndex}, familyName::S, tenor::TenorPeriod, fixingDays::I, currency::AbstractCurrency, fixingCalendar::B,
  #                               convention::C, endOfMonth::Bool, dc::DC) =
  #   new{S, I, B, C, DC, TermStructure}(familyName, tenor, fixingDays, currency, fixingCalendar, convention, endOfMonth, dc)
  #
  # call{S, I, B, C, DC, T}(::Type{IborIndex}, familyName::S, tenor::TenorPeriod, fixingDays::I, currency::AbstractCurrency, fixingCalendar::B,
  #                               convention::C, endOfMonth::Bool, dc::DC, ts::T) =
  #   new{S, I, B, C, DC, T}(familyName, tenor, fixingDays, currency, fixingCalendar, convention, endOfMonth, dc, ts)
  # IborIndex(familyName::S, tenor::TenorPeriod, fixingDays::I, currency::AbstractCurrency, fixingCalendar::B, convention::C,
  #                   endOfMonth::Bool, dc::DC, ts::T = NullTermStructure()) = IborIndex(familyName, tenor, fixingDays, currency, fixingCalendar, convention,
  #                                                                                           endOfMonth, dc, ts)
end

IborIndex{TP <: TenorPeriod, CUR <: AbstractCurrency, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure}(familyName::String, tenor::TP, fixingDays::Int, currency::CUR,
          fixingCalendar::B, convention::C, endOfMonth::Bool, dc::DC, ts::T = NullTermStructure()) =
          IborIndex{TP, CUR, B, C, DC, T}(familyName, tenor, fixingDays, currency, fixingCalendar, convention, endOfMonth, dc, ts, Dict{Date, Float64}())


immutable LiborIndex{TP <: TenorPeriod, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure} <: InterestRateIndex
  familyName::String
  tenor::TP
  fixingDays::Int
  currency::Currency
  fixingCalendar::B
  jointCalendar::JointCalendar
  convention::C
  endOfMonth::Bool
  dc::DC
  ts::T
  pastFixings::Dict{Date, Float64}

  # call{S, I, B, C, DC}(::Type{LiborIndex}, familyName::S, tenor::TenorPeriod, fixingDays::I, currency::Currency, fixingCalendar::B,
  #                               jointCalendar::JointCalendar, convention::C, endOfMonth::Bool, dc::DC) =
  #   new{S, I, B, C, DC, TermStructure}(familyName, tenor, fixingDays, currency, fixingCalendar, jointCalendar, convention, endOfMonth, dc)

end

# catch all constructor #
LiborIndex{TP <: TenorPeriod, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount, T <: TermStructure}(familyName::String, tenor::TP, fixingDays::Int, currency::Currency, fixingCalendar::B,
                  jointCalendar::JointCalendar, convention::C, endOfMonth::Bool, dc::DC, ts::T = NullTermStructure()) = LiborIndex{TP, B, C, DC, T}(familyName, tenor, fixingDays, currency, fixingCalendar, jointCalendar, convention,
                                                                                          endOfMonth, dc, ts, Dict{Date, Float64}())


function LiborIndex(familyName::String, tenor::TenorPeriod, fixingDays::Int, currency::Currency, fixingCalendar::BusinessCalendar, dc::DayCount, yts::YieldTermStructure)
  endOfMonth = libor_eom(tenor.period)
  conv = libor_conv(tenor.period)
  jc = JointCalendar(QuantLib.Time.UKLSECalendar(), fixingCalendar)

  return LiborIndex(familyName, tenor, fixingDays, currency, fixingCalendar, jc, conv, endOfMonth, dc, yts)
end



fixing_date(idx::InterestRateIndex, d::Date) = advance(Dates.Day(-idx.fixingDays), idx.fixingCalendar, d, idx.convention)
maturity_date(idx::IborIndex, d::Date) = advance(idx.tenor.period, idx.fixingCalendar, d, idx.convention)
value_date(idx::InterestRateIndex, d::Date) = advance(Dates.Day(idx.fixingDays), idx.fixingCalendar, d)

function fixing(idx::InterestRateIndex, ts::TermStructure, _fixing_date::Date, forecast_todays_fixing::Bool=true)
  today = settings.evaluation_date
  if _fixing_date > today || (_fixing_date == today && forecast_todays_fixing)
    return forecast_fixing(idx, ts, _fixing_date)
  end

  # this is a path fixing
  pastFix = get(idx.pastFixings, _fixing_date, -1.0)

  if pastFix == -1.0
    # none found, default to forecast
    return forecast_fixing(idx, ts, _fixing_date)
  else
    return pastFix
  end
  # error("Not yet implemented for older dates than eval date")
end

function forecast_fixing(idx::InterestRateIndex, ts::TermStructure, _fixing_date::Date)
  d1 = value_date(idx, _fixing_date)
  d2 = maturity_date(idx, d1)
  t = year_fraction(idx.dc, d1, d2)
  return forecast_fixing(idx, ts, d1, d2, t)
end

function forecast_fixing(idx::InterestRateIndex, ts::TermStructure, d1::Date, d2::Date, t::Float64)
  disc1 = discount(ts, d1)
  disc2 = discount(ts, d2)

  return (disc1 / disc2 - 1.0) / t
end

is_valid_fixing_date(idx::InterestRateIndex, d::Date) = is_business_day(idx.fixingCalendar, d)

function add_fixing!(idx::InterestRateIndex, d::Date, fixingVal::Float64)
  idx.pastFixings[d] = fixingVal

  return idx
end

# Libor methods
function value_date(idx::LiborIndex, d::Date)
  new_d = advance(Dates.Day(idx.fixingDays), idx.fixingCalendar, d, idx.convention)
  return adjust(idx.jointCalendar, idx.convention, new_d)
end

maturity_date(idx::LiborIndex, d::Date) = advance(idx.tenor.period, idx.jointCalendar, d, idx.convention)

# types of indexes
euribor_index(tenor::TenorPeriod) = IborIndex("Euribor", tenor, 2, EURCurrency(), QuantLib.Time.TargetCalendar(), euribor_conv(tenor.period), euribor_eom(tenor.period), QuantLib.Time.Actual360())
euribor_index(tenor::TenorPeriod, ts::TermStructure) = IborIndex("Euribor", tenor, 2, EURCurrency(), QuantLib.Time.TargetCalendar(), euribor_conv(tenor.period), euribor_eom(tenor.period), QuantLib.Time.Actual360(), ts)

function usd_libor_index(tenor::TenorPeriod, yts::YieldTermStructure)
  return LiborIndex("USDLibor", tenor, 2, USDCurrency(), QuantLib.Time.USSettlementCalendar(), QuantLib.Time.Actual360(), yts)
end

euribor_conv(::Union{Base.Dates.Day, Base.Dates.Week}) = QuantLib.Time.Following()
euribor_conv(::Union{Base.Dates.Month, Base.Dates.Year}) = QuantLib.Time.ModifiedFollowing()

euribor_eom(::Union{Base.Dates.Day, Base.Dates.Week}) = false
euribor_eom(::Union{Base.Dates.Month, Base.Dates.Year}) = true

libor_conv(::Union{Base.Dates.Day, Base.Dates.Week}) = QuantLib.Time.Following()
libor_conv(::Union{Base.Dates.Month, Base.Dates.Year}) = QuantLib.Time.ModifiedFollowing()

libor_eom(::Union{Base.Dates.Day, Base.Dates.Week}) = false
libor_eom(::Union{Base.Dates.Month, Base.Dates.Year}) = true

# clone methods #
clone{TP, CUR, B, C, DC, T}(idx::IborIndex{TP, CUR, B, C, DC, T}, ts::TermStructure = idx.ts) =
      IborIndex(idx.familyName, idx.tenor, idx.fixingDays, idx.currency, idx.fixingCalendar, idx.convention, idx.endOfMonth, idx.dc, ts)
clone(idx::LiborIndex, ts::TermStructure = idx.ts) = LiborIndex(idx.familyName, idx.tenor, idx.fixingDays, idx.currency, idx.fixingCalendar, idx.jointCalendar, idx.convention, idx.endOfMonth, idx.dc, ts)

update_termstructure(idx::InterestRateIndex, ts::TermStructure) = clone(idx, ts)
