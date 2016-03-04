# Null
type FdmNullStepCondition <: StepCondition end

apply_to!(cond::FdmNullStepCondition, a::Vector{Float64}, ::Float64) = cond, a

type FdmDividendHandler{F <: FdmMesher, I <: Integer} <: StepCondition
  x::Vector{Float64}
  dividendTimes::Vector{Float64}
  dividendDates::Vector{Date}
  dividends::Vector{Float64}
  mesher::F
  equityDirection::I
end

function FdmDividendHandler(schedule::DividendSchedule, mesher::FdmMesher, refDate::Date, dc::DayCount, equityDirection::Int)
  schedLength = length(schedule.dividends)
  dividends = zeros(schedLength)
  dividendDates = Vector{Date}(schedLength)
  dividendTimes = zeros(schedLength)
  x = zeros(mesher.layout.dim[equityDirection])

  for i = 1:schedLength
    dividends = amount(schedule.dividends[i])
    dividendDates = date(schedule.dividends[i])
    dividendTimes = year_fraction(dc, refDate, date(schedule.dividends[i]))
  end

  tmp = get_locations(mesher, equityDirection)
  spacing = mesher.layout.spacing[equityDirection]

  for i = 1:length(x)
    iter = ((i - 1) * spacing) + 1
    x[i] = exp(tmp[iter])
  end

  return FdmDividendHandler(x, dividendTimes, dividendDates, dividends, mesher, equityDirection)
end

type FdmAmericanStepCondition{F <: FdmMesher, C <: FdmInnerValueCalculator} <: StepCondition
  mesher::F
  calculator::C
end

type FdmBermudanStepCondition{F <: FdmMesher, C <: FdmInnerValueCalculator} <: StepCondition
  mesher::F
  calculator::C
  exerciseTimes::Vector{Float64}
end

function FdmBermudanStepCondition(exerciseDates::Vector{Date}, refDate::Date, dc::DayCount, mesher::FdmMesher, calculator::FdmInnerValueCalculator)
  exerciseTimes = zeros(length(exerciseDates))
  for i = 1:length(exerciseTimes)
    exerciseTimes[i] = year_fraction(dc, refDate, exerciseDates[i])
  end

  return FdmBermudanStepCondition(mesher, calculator, exerciseTimes)
end

function apply_to!(cond::FdmBermudanStepCondition, a::Vector{Float64}, t::Float64)
  if findfirst(cond.exerciseTimes, t) != 0
    layout = cond.mesher.layout
    dims = length(layout.dim)
    coords = ones(Int, dims)

    locations = zeros(dims)

    for i = 1:layout.size
      for j = 1:dims
        locations[dims] = get_location(cond.mesher, coords, j)
      end

      innerValue, cond.calculator = inner_value(cond.calculator, coords, i, t)
      if innerValue > a[i]
        a[i] = innerValue
      end

      iter_coords!(coords, cond.mesher.layout.dim)
    end
  end

  return cond, a
end

type FdmSnapshotCondition <: StepCondition
  t::Float64
  a::Vector{Float64}
end

FdmSnapshotCondition(t::Float64) = FdmSnapshotCondition(t, Vector{Float64}(0))

function apply_to!(cond::FdmSnapshotCondition, a::Vector{Float64}, t::Float64)
  if cond.t == t
    cond.a = a
  end

  return cond, a
end

type FdmStepConditionComposite{C <: StepCondition} <: StepCondition
  stoppingTimes::Vector{Float64}
  conditions::Vector{C}
end

# Constructors
function vanilla_FdmStepConditionComposite(cashFlow::DividendSchedule, exercise::Exercise, mesher::FdmMesher, calculator::FdmInnerValueCalculator,
                                          refDate::Date, dc::DayCount)
  stoppingTimes = Vector{Float64}()
  stepConditions = Vector{StepCondition}()

  if length(cashFlow.dividends) > 0
    dividendCondition = FdmDividendHandler(cashFlow, mesher, refDate, dc, 1)
    push!(stepConditions, dividendCondition)
    append!(stoppingTimes, dividendCondition.dividendTimes)
  end

  if isa(exercise, AmericanExercise)
    push!(stepConditions, FdmAmericanStepCondition(mesher, calculator))
  elseif isa(exercise, BermudanExercise)
    bermudanCondition = FdmBermudanStepCondition(exercise.dates, refDate, dc, mesher, calculator)
    push!(stepConditions, bermudanCondition)
    append!(stoppingTimes, bermudanCondition.exerciseTimes)
  end

  return FdmStepConditionComposite(stoppingTimes, stepConditions)
end

function join_conditions_FdmStepConditionComposite(c1::FdmSnapshotCondition, c2::FdmStepConditionComposite)
  stoppingTimes = c2.stoppingTimes
  push!(stoppingTimes, c1.t)

  conditions = Vector{StepCondition}(2)
  conditions[1] = c2
  conditions[2] = c1

  return FdmStepConditionComposite(stoppingTimes, conditions)
end

function apply_to!(cond::FdmStepConditionComposite, a::Vector{Float64}, t::Float64)
  for c in cond.conditions
    apply_to!(c, a, t)
  end

  return cond, a
end

type ArrayWrapper <: CurveWrapper
  value::Vector{Float64}
end

get_value(wrapper::ArrayWrapper, ::Vector{Float64}, idx::Int) = wrapper.value[idx]

type PayoffWrapper{P <: StrikedTypePayoff} <: CurveWrapper
  payoff::P
end

type AmericanStepConditionType <: StepType end

type CurveDependentStepCondition{T <: StepType, C <: CurveWrapper} <: StepCondition
  stepType::T
  curveItem::C
end

typealias AmericanStepCondition{C} CurveDependentStepCondition{AmericanStepConditionType, C}

build_AmericanStepCondition(intrinsicValues::Vector{Float64}) = CurveDependentStepCondition(AmericanStepConditionType(), ArrayWrapper(intrinsicValues))

apply_to_value(::AmericanStepCondition, current::Float64, intrinsic::Float64) = max(current, intrinsic)

function apply_to!(cond::CurveDependentStepCondition, a::Vector{Float64}, ::Float64)
  for i in eachindex(a)
    a[i] = apply_to_value(cond, a[i], get_value(cond, a, i))
  end

  return cond, a
end

get_value(cond::CurveDependentStepCondition, a::Vector{Float64}, idx::Int) = get_value(cond.curveItem, a, idx)

type StepConditionSet{C <: StepCondition} <: StepCondition
  conditions::Vector{C}
end

function apply_to!(cond::StepConditionSet, a::Vector{Vector{Float64}}, t::Float64)
  for i in eachindex(cond.conditions)
    apply_to!(cond.conditions[i], a[i], t)
  end

  return cond, a
end
