struct AmericanExercise <: EarlyExercise
  dates::Vector{Date}
end

struct BermudanExercise <: EarlyExercise
  dates::Vector{Date}
end

struct EuropeanExercise <: Exercise
  dates::Vector{Date}
end

EuropeanExercise(d::Date) = EuropeanExercise([d])
AmericanExercise(d1::Date, d2::Date) = AmericanExercise(Date[d1, d2])

struct RebatedExercise{E <: Exercise, C <: BusinessCalendar, B <: BusinessDayConvention} <: Exercise
  exercise::E
  rebate::Float64
  rebateSettlementDays::Int
  rebatePaymentCalendar::C
  rebatePaymentConvention::B

  RebatedExercise{E, C, B}(exercise::E,
                          rebate::Float64 = 0.0,
                          rebateSettlementDays::Int = 0,
                          rebatePaymentCalendar::C = NullCalendar(),
                          rebatePaymentConvention::B = Following()) where {E, C, B} =
                new{E, C, B}(exercise, rebate, rebateSettlementDays, rebatePaymentCalendar, rebatePaymentConvention)
end
