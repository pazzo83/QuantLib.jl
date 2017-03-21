type FdHullWhiteSwaptionEngine{Y <: TermStructure, F <: FdmSchemeDescType} <: PricingEngine
  model::HullWhite{AffineModelType, Y}
  tGrid::Int
  xGrid::Int
  dampingSteps::Int
  invEps::Float64
  schemeDesc::FdmSchemeDesc{F}
  ts::Y
end

FdHullWhiteSwaptionEngine{F <: FdmSchemeDescType}(model::HullWhite, tGrid::Int = 100, xGrid::Int = 100, dampingSteps::Int = 0, invEps::Float64 = 1e-5,
                  schemeDesc::FdmSchemeDesc{F} = FdmSchemeDesc(Douglas())) =
                  FdHullWhiteSwaptionEngine{typeof(model.ts), F}(model, tGrid, xGrid, dampingSteps, invEps, schemeDesc, model.ts)
# methods #
function _calculate!(pe::FdHullWhiteSwaptionEngine, swaption::Swaption)
  # 1. Term structure
  ts = pe.ts

  # 2. Mesher
  dc = ts.dc
  refDate = reference_date(ts)
  maturity = year_fraction(dc, refDate, swaption.exercise.dates[end])

  process = OrnsteinUhlenbeckProcess(get_a(pe.model), get_sigma(pe.model))
  shortRateMesher = FdmSimpleProcess1dMesher(pe.xGrid, process, maturity, 1, pe.invEps)
  mesher = FdmMesherComposite(shortRateMesher)

  # 3. Inner Value calculator
  exerciseDates = swaption.exercise.dates
  t2d = Dict{Float64, Date}()
  @simd for i = 1:length(exerciseDates)
    @inbounds t = year_fraction(dc, refDate, exerciseDates[i])
    @inbounds t2d[t] = exerciseDates[i]
  end

  disTs = pe.model.ts
  fwdTs = swaption.swap.iborIndex.ts

  # TODO check that day counts and ref dates match btwn the two term structures
  fwdModel = HullWhite(fwdTs, get_a(pe.model), get_sigma(pe.model))
  calculator = FdmAffineModelSwapInnerValue(pe.model, fwdModel, swaption.swap, t2d, mesher, 1)

  # 4. Step Conditions
  conditions = vanilla_FdmStepConditionComposite(DividendSchedule(), swaption.exercise, mesher, calculator, refDate, dc)

  # 5. Boundary conditions
  boundaries = FdmBoundaryConditionSet()

  # 6. Solver
  solverDesc = FdmSolverDesc(mesher, boundaries, conditions, calculator, maturity, pe.tGrid, pe.dampingSteps)

  solver = FdmHullWhiteSolver(pe.model, solverDesc, pe.schemeDesc)
  swaption.results.value = value_at(solver, 0.0)
end
