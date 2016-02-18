type FdG2SwaptionEngine{I <: Integer, Y <: TermStructure} <: PricingEngine
  model::G2
  tGrid::I
  xGrid::I
  yGrid::I
  dampingSteps::I
  invEps::Float64
  schemeDesc::FdmSchemeDesc
  ts::Y
end

FdG2SwaptionEngine(model::G2, tGrid::Int = 100, xGrid::Int = 50, yGrid::Int = 50, dampingSteps::Int = 0, invEps::Float64 = 1e-5,
                  schemeDesc::FdmSchemeDesc = FdmSchemeDesc(Hundsdorfer())) =
                  FdG2SwaptionEngine(model, tGrid, xGrid, yGrid, dampingSteps, invEps, schemeDesc, model.ts)

# methods #
function _calculate!(pe::FdG2SwaptionEngine, swaption::Swaption)
  # 1. Term structure
  ts = pe.ts

  # 2. Mesher
  dc = ts.dc
  refDate = reference_date(ts)
  maturity = year_fraction(dc, refDate, swaption.exercise.dates[end])

  process1 = OrnsteinUhlenbeckProcess(get_a(pe.model), get_sigma(pe.model))
  process2 = OrnsteinUhlenbeckProcess(get_b(pe.model), get_eta(pe.model))

  xMesher = FdmSimpleProcess1dMesher(pe.xGrid, process1, maturity, 1, pe.invEps)
  yMesher = FdmSimpleProcess1dMesher(pe.yGrid, process2, maturity, 1, pe.invEps)
  mesher = FdmMesherComposite(xMesher, yMesher)

  # 3. Inner Value calculator
  exerciseDates = swaption.exercise.dates
  t2d = Dict{Float64, Date}()
  for i = 1:length(exerciseDates)
    t = year_fraction(dc, refDate, exerciseDates[i])
    t2d[t] = exerciseDates[i]
  end

  disTs = pe.model.ts
  fwdTs = swaption.swap.iborIndex.ts

  fwdModel = G2(fwdTs, get_a(pe.model), get_sigma(pe.model), get_b(pe.model), get_eta(pe.model), get_rho(pe.model))
  calculator = FdmAffineModelSwapInnerValue(pe.model, fwdModel, swaption.swap, t2d, mesher, 1)

  # 4. Step Conditions
  conditions = vanilla_FdmStepConditionComposite(DividendSchedule(), swaption.exercise, mesher, calculator, refDate, dc)

  # 5. Boundary conditions
  boundaries = FdmBoundaryConditionSet()

  # 6. Solver
  solverDesc = FdmSolverDesc(mesher, boundaries, conditions, calculator, maturity, pe.tGrid, pe.dampingSteps)
  solver = FdmG2Solver(pe.model, solverDesc, pe.schemeDesc)
  swaption.results.value = value_at(solver, 0.0, 0.0)
end
