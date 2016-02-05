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
