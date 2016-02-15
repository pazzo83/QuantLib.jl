type AmericanExercise <: EarlyExercise
  dates::Vector{Date}
end

type BermudanExercise <: EarlyExercise
  dates::Vector{Date}
end

type EuropeanExercise <: Exercise
  dates::Vector{Date}
end

EuropeanExercise(d::Date) = EuropeanExercise([d])

type RebatedExercise{E <: Exercise, I <: Integer, C <: BusinessCalendar, B <: BusinessDayConvention} <: Exercise
  exercise::E
  rebate::Float64
  rebateSettlementDays::I
  rebatePaymentCalendar::C
  rebatePaymentConvention::B

  RebatedExercise(exercise::E, rebate::Float64 = 0.0, rebateSettlementDays::I = 0, rebatePaymentCalendar::C = NullCalendar(),
                rebatePaymentConvention::B = Following()) = new{E, I, C, B}(exercise, rebate, rebateSettlementDays, rebatePaymentCalendar, rebatePaymentConvention)
end
