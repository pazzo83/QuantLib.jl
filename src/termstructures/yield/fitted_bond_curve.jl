type FittedBondDiscountCurve{C <: BusinessCalendar, B <: BondHelper, DC <: DayCount, F <: FittingMethod} <: Curve
  lazyMixin::LazyMixin
  settlementDays::Int
  referenceDate::Date
  calendar::C
  bonds::Vector{B}
  dc::DC
  fittingMethod::F
  accuracy::Float64
  maxEvaluations::Int
  simplexLambda::Float64

  FittedBondDiscountCurve(settlementDays::Int,
                          referenceDate::Date,
                          calendar::C,
                          bonds::Vector{B},
                          dc::DC,
                          fittingMethod::F,
                          accuracy::Float64,
                          maxEvaluations::Int,
                          simplexLambda::Float64) =

                          (x = new(LazyMixin(), settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda);
                          x.fittingMethod.commons.costFunction.curve = x)

  # FittedBondDiscountCurve(settlementDays::Integer, referenceDate::Date, calendar::BusinessCalendar, bonds::Vector{B}, dc::DayCount, fittingMethod::FittingMethod, accuracy::Float64=1e-10,
  #                                     maxEvaluations::Integer=10000, simplexLambda::Float64=1.0) =
  #                         new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)
    #n = length(bonds)
  #   println("hi")
  #   this = new(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)
  #   println("also hi")
  #   this.fittingMethod.costFunction.curve = this
  #
  #   return this
  # end
end

FittedBondDiscountCurve{C <: BusinessCalendar, B <: BondHelper, DC <: DayCount, F <: FittingMethod}(settlementDays::Int, referenceDate::Date, calendar::C, bonds::Vector{B}, dc::DC, fittingMethod::F, accuracy::Float64=1e-10,
                                     maxEvaluations::Int=10000, simplexLambda::Float64=1.0) =
                                     FittedBondDiscountCurve{C, B, DC, F}(settlementDays, referenceDate, calendar, bonds, dc, fittingMethod, accuracy, maxEvaluations, simplexLambda)

### Fitted curve methods ###
discount_impl(curve::FittedBondDiscountCurve, t::Float64) = discount_function(curve.fittingMethod, curve.fittingMethod.commons.solution, t)

function initialize!(curve::FittedBondDiscountCurve)
 # yield conventions
 dc = curve.dc
 yield_comp = CompoundedCompounding()
 freq = Annual()

 n = length(curve.bonds)
 cost_f = curve.fittingMethod.commons.costFunction

 squared_sum = 0.0
 @inbounds for i in eachindex(curve.bonds)
   bond = curve.bonds[i].bond
   leg = bond.cashflows
   clean_price = bond.faceAmount
   bond_settlement = get_settlement_date(bond)

   # get the ytm of the bond
   ytm = yield(bond, clean_price, dc, yield_comp, freq, bond_settlement)
   dur = duration(bond, ytm, dc, yield_comp, freq, ModifiedDuration(), bond_settlement)

   curve.fittingMethod.commons.weights[i] = 1.0 / dur
   squared_sum += curve.fittingMethod.commons.weights[i] * curve.fittingMethod.commons.weights[i]

   cf = bond.cashflows
   for k = 1:length(cf.coupons) # for redemption
     if !has_occurred(cf.coupons[k], bond_settlement)
       cost_f.firstCashFlow[i] = k
       break
     end
   end
 end

 divide_array_by_self!(curve.fittingMethod.commons.weights, sqrt(squared_sum))

 return curve
end

function perform_calculations!(curve::FittedBondDiscountCurve)
 cost_f = curve.fittingMethod.commons.costFunction
 constraint = NoConstraint()

 x = zeros(curve.fittingMethod.size)
 # x = Vector{DD}(curve.fittingMethod.size)

 if length(curve.fittingMethod.commons.guessSolution) > 0
   x = curve.fittingMethod.commons.guessSolution
 end

 simplex = Simplex(curve.simplexLambda)
 problem = Problem(cost_f, constraint, x)

 max_stationary_state_iterations = 100
 root_epsilon = curve.accuracy
 function_epsilon = curve.accuracy
 gradient_norm_epsilon = curve.accuracy

 end_criteria = EndCriteria(curve.maxEvaluations, max_stationary_state_iterations, root_epsilon, function_epsilon, gradient_norm_epsilon)

 minimize!(simplex, problem, end_criteria)
 curve.fittingMethod.commons.solution = problem.currentValue

 number_of_iterations = problem.functionEvaluation
 cost_value = problem.functionValue

 curve.fittingMethod.commons.guessSolution = curve.fittingMethod.commons.solution
 curve.fittingMethod.commons.numberOfIterations = number_of_iterations
 curve.fittingMethod.commons.minimumCostValue = cost_value

 # we have calculated

 return curve
end
