type LongstaffSchwartzExerciseStrategy{M <: MarketModelBasisSystem, E <: MarketModelExerciseValue, E2 <: MarketModelExerciseValue} <: ExerciseStrategy
  basisSystem::M
  basisCoefficients::Vector{Vector{Float64}}
  exercise::E
  control::E2
  numeraires::Vector{Int}
  currentIndex::Int
  principalInNumerairePortfolio::Float64
  newPrincipal::Float64
  exerciseTimes::Vector{Float64}
  relevantTimes::Vector{Float64}
  isBasisTime::BitArray{1}
  isRebateTime::BitArray{1}
  isControlTime::BitArray{1}
  isExerciseTime::BitArray{1}
  rebateDiscounters::Vector{MarketModelDiscounter}
  controlDiscounters::Vector{MarketModelDiscounter}
  basisValues::Vector{Vector{Float64}}
  exerciseIndex::Vector{Int}
end

function LongstaffSchwartzExerciseStrategy(basisSystem::MarketModelBasisSystem,
                                          basisCoefficients::Vector{Vector{Float64}},
                                          evolution::EvolutionDescription,
                                          numeraires::Vector{Int},
                                          exercise::MarketModelExerciseValue,
                                          control::MarketModelExerciseValue)

  basisSystem = clone(basisSystem)
  exercise = clone(exercise)
  control = clone(control)

  check_compatibility(evolution, numeraires)
  relevantTimes = evolution.evolutionTimes
  isBasisTime = is_in_subset(relevantTimes, basisSystem.evolution.evolutionTimes)
  isRebateTime = is_in_subset(relevantTimes, exercise.evolution.evolutionTimes)
  isControlTime = is_in_subset(relevantTimes, control.evolution.evolutionTimes)

  exerciseIndex = Vector{Int}(length(relevantTimes))
  isExerciseTime = falses(length(exerciseIndex))
  exerciseTimes = Vector{Float64}()
  v = exercise.isExerciseTime
  exercises = idx = 1
  for i in eachindex(relevantTimes)
    exerciseIndex[i] = exercises
    if isRebateTime[i]
      isExerciseTime[i] = v[idx]
      idx += 1
      if isExerciseTime[i]
        push!(exerciseTimes, relevantTimes[i])
        exercises += 1
      end
    end
  end

  rateTimes = evolution.rateTimes
  rebateTimes = possible_cash_flow_times(exercise)
  rebateDiscounters = MarketModelDiscounter[MarketModelDiscounter(rebateTimes[i], rateTimes) for i in eachindex(rebateTimes)]

  controlTimes = possible_cash_flow_times(control)
  controlDiscounters = MarketModelDiscounter[MarketModelDiscounter(controlTimes[i], rateTimes) for i in eachindex(controlTimes)]

  basisSizes = number_of_functions(basisSystem)
  basisValues = Vector{Float64}[Vector{Float64}(basisSizes[i]) for i in eachindex(basisSizes)]

  return LongstaffSchwartzExerciseStrategy(basisSystem, basisCoefficients, exercise, control, numeraires, -1, 0.0, 0.0, exerciseTimes, relevantTimes,
                                          isBasisTime, isRebateTime, isControlTime, isExerciseTime, rebateDiscounters, controlDiscounters,
                                          basisValues, exerciseIndex)
end

relevant_times(ls::LongstaffSchwartzExerciseStrategy) = ls.relevantTimes

function reset!(ls::LongstaffSchwartzExerciseStrategy)
  reset!(ls.exercise)
  reset!(ls.control)
  reset!(ls.basisSystem)
  ls.currentIndex = 1
  ls.principalInNumerairePortfolio = ls.newPrincipal = 1.0

  return ls
end

function next_step!(ls::LongstaffSchwartzExerciseStrategy, currentState::CurveState)
  ls.principalInNumerairePortfolio = ls.newPrincipal

  if ls.isRebateTime[ls.currentIndex]
    next_step!(ls.exercise, currentState)
  end

  if ls.isControlTime[ls.currentIndex]
    next_step!(ls.control, currentState)
  end

  if ls.isBasisTime[ls.currentIndex]
    next_step!(ls.basisSystem, currentState)
  end

  if ls.currentIndex < length(ls.numeraires)
    numeraire = ls.numeraires[ls.currentIndex]
    nextNumeraire = ls.numeraires[ls.currentIndex + 1]
    ls.newPrincipal *= discount_ratio(currentState, numeraire, nextNumeraire)
  end

  ls.currentIndex += 1

  return ls
end

function get_exercise!(ls::LongstaffSchwartzExerciseStrategy, currentState::CurveState)
  exerciseIndex = ls.exerciseIndex[ls.currentIndex-1]
  # println("IDX: ", exerciseIndex)

  exerciseCF = get_value(ls.exercise, currentState)
  exerciseValue = exerciseCF.amount * numeraire_bonds(ls.rebateDiscounters[exerciseCF.timeIndex], currentState, ls.numeraires[ls.currentIndex - 1]) /
                  ls.principalInNumerairePortfolio

  controlCF = get_value(ls.control, currentState)
  controlValue = controlCF.amount * numeraire_bonds(ls.controlDiscounters[controlCF.timeIndex], currentState, ls.numeraires[ls.currentIndex - 1]) /
                  ls.principalInNumerairePortfolio

  set_values!(ls.basisSystem, currentState, ls.basisValues[exerciseIndex])

  alphas = ls.basisCoefficients[exerciseIndex]
  continuationValue = dot(alphas, ls.basisValues[exerciseIndex])

  # println("alphas: ", alphas)
  # println("basis: ", ls.basisValues[exerciseIndex])
  # println("exerciseVal: ", exerciseValue)
  # println("continuationValue: ", continuationValue)

  return exerciseValue >= continuationValue
end

clone(ls::LongstaffSchwartzExerciseStrategy) = LongstaffSchwartzExerciseStrategy(clone(ls.basisSystem), deepcopy(ls.basisCoefficients), clone(ls.exercise),
                                                clone(ls.control), copy(ls.numeraires), ls.currentIndex, ls.principalInNumerairePortfolio, ls.newPrincipal,
                                                copy(ls.exerciseTimes), copy(ls.relevantTimes), copy(ls.isBasisTime), copy(ls.isRebateTime), copy(ls.isControlTime),
                                                copy(ls.isExerciseTime), deepcopy(ls.rebateDiscounters), deepcopy(ls.controlDiscounters), deepcopy(ls.basisValues),
                                                copy(ls.exerciseIndex))
