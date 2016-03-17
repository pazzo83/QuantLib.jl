using QuantLib

function InverseFloater(rateLevel::Float64)
  numberOfRates = 20
  accrual = 0.5
  firstTime = 0.5

  strike = 0.15
  fixedMultiplier = 2.0
  floatingSpread = 0.0
  payer = true

  rateTimes = Float64[firstTime + (i-1) * accrual for i = 1:numberOfRates + 1]

  accruals = fill(accrual, numberOfRates)
  fixedStrikes = fill(strike, numberOfRates)
  floatingSpreads = fill(floatingSpread, numberOfRates)
  fixedMultipliers = fill(fixedMultiplier, numberOfRates)

  paymentTimes = Float64[firstTime + i * accrual for i = 1:numberOfRates]

  inverseFloater = MultiStepInverseFloater(rateTimes, accruals, accruals, fixedStrikes, fixedMultipliers, floatingSpreads, paymentTimes, payer)

  # exercise schedule, we can exercise on any rate time except the last one
  exerciseTimes = rateTimes[1:end-1]

  # naive exercise strategy, exercise above a trigger level
  trigger = 0.05
  swapTriggers = fill(trigger, length(exerciseTimes))
  naifStrategy = SwapRateTrigger(rateTimes, swapTriggers, exerciseTimes)

  # Longstaff-Schwartz exercise strategy
  collectedData = Vector{NodeData}()
  basisCoefficients = Vector{Float64}()

  # control that does nothing, need it because some control is required
  control = NothingExerciseValue(rateTimes)
end

function main()
  InverseFloater(5 / 100.0)
end
