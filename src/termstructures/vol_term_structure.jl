type NullOptionVolatilityStructure <: OptionletVolatilityStructure end

type ConstantOptionVolatility{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount} <: OptionletVolatilityStructure
  settlementDays::I
  referenceDate::Date
  calendar::B
  bdc::C
  volatility::Float64
  dc::DC
end

function ConstantOptionVolatility{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount}(settlementDays::I, calendar::B, bdc::C, volatility::Float64, dc::DC)
  today = settings.evaluation_date
  ref_date = advance(Dates.Day(settlementDays), calendar, today, bdc)
  ConstantOptionVolatility(settlementDays, ref_date, calendar, bdc, volatility, dc)
end

# Swaption Volatility structures
type ConstantSwaptionVolatility{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount} <: SwaptionVolatilityStructure
  settlementDays::I
  referenceDate::Date
  calendar::B
  bdc::C
  volatility::Quote
  dc::DC
end

# floating reference date, floating market data
ConstantSwaptionVolatility{I <: Integer, B <: BusinessCalendar, C <: BusinessDayConvention, DC <: DayCount}(settlementDays::I, cal::B, bdc::C, volatility::Quote, dc::DC) =
                          ConstantSwaptionVolatility(settlementDays, Date(), cal, bdc, volatility, dc)

function black_varience(ovs::OptionletVolatilityStructure, option_date::Date, strike::Float64)
  v = calc_volatility(ovs, option_date, strike)
  t = time_from_reference(ovs, option_date)
  return v * v * t
end

function calc_volatility(ovs::OptionletVolatilityStructure, option_date::Date, strike::Float64)
  # TODO checks - see optionletvolatilitystructure.hpp, volatility
  return volatility_impl(ovs, option_date, strike)
end

volatility_impl(ovs::OptionletVolatilityStructure, option_date::Date, strike::Float64) =
              volatility_impl(ovs, time_from_reference(ovs, option_date), strike)

volatility_impl(const_opt_vol::ConstantOptionVolatility, ::Float64, ::Float64) = const_opt_vol.volatility

# Swaption Volatility methods
function swap_length{V <: SwaptionVolatilityStructure}(swapVol::V, _start::Date, _end::Date)
  result = round(float(_end - _start) / 365.25 * 12.0) # month unit
  result /= 12.0 # year unit
  return result
end

function calc_volatility{V <: SwaptionVolatilityStructure}(swapVol::V, optionDate::Date, swapLength::Float64, strike::Float64)
  # TODO add checks
  optionTime = time_from_reference(swapVol, optionDate)
  return volatility_impl(swapVol, optionTime, swapLength, strike)
end

volatility_impl(swapVol::ConstantSwaptionVolatility, ::Float64, ::Float64, ::Float64) = swapVol.volatility.value

function black_varience{V <: SwaptionVolatilityStructure}(swapVol::V, optionDate::Date, swapLength::Float64, strike::Float64)
  v = calc_volatility(swapVol, optionDate, swapLength, strike)
  optionTime = time_from_reference(swapVol, optionDate)
  return v * v * optionTime
end
