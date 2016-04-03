# QuantLib
module QuantLib

# functions overridden from based
import Base.findprev, Base.findnext

function findprev(testf::Function, A, start::Integer, val)
  for i = start:-1:1
    testf(A[i], val) && return i
  end
  0
end

function findnext(testf::Function, A, start::Integer, val)
  for i = start:length(A)
    if testf(A[i], val)
      return i
    end
  end
  return 0
end

function upper_bound{T}(vec::Vector{T}, x::T)
  if x > vec[end]
    return length(vec) + 1
  end
  found = searchsortedlast(vec, x)
  if found != length(vec)
    found += 1
  end

  return found
end

# Math module
include("math/Math.jl")

# Time module
include("time/Time.jl")

# MAIN MODULE CODE
using QuantLib.Math, QuantLib.Time
using Distributions

# Constants #
const basisPoint = 0.0001
const JQ_MAX = 1.79769e308
const JQ_MIN = 2.22507e-308
const M_SQRT2 = 1.41421356237309504880168872420969808   # sqrt(2)
const M_SQRT1_2 = 0.707106781186547524400844362104849039  # 1/sqrt(2)
const M_SQRTPI = 1.77245385090551602792981
const M_SQRT_2 = 0.7071067811865475244008443621048490392848359376887
const M_1_SQRTPI = 0.564189583547756286948

export
    # abstract_types.jl
    LazyObject,

    Exercise, EarlyExercise, CompoundingType, TermStructure, YieldTermStructure, InterpolatedCurve, BootstrapTrait, Bootstrap, BootstrapHelper, BondHelper, RateHelper,
    FittingMethod, CashFlows, CashFlow, Coupon, CouponPricer, IborCouponPricer, Instrument, Bond, Swap, SwapType, PricingEngine, Duration, AbstractRate, Results,
    InterestRateIndex, AbstractCurrency, Parameter, CalibrationHelper, ShortRateModel, FdmMesher, OptionType, AbstractMarketModel, OptionletVolatilityStructure,

    # lazy.jl
    LazyMixin, calculate!, recalculate!,

    # Methods
    # methods/lattice.jl
    TreeLattice1D, Branching, TrinomialTree,

    # methods/binomial_tree.jl
    JarrowRudd, CoxRossRubinstein, AdditiveEQP, Trigeorgis, Tian, LeisenReimer, Joshi4,

    # methods/monte_carlo/utilities.jl
    NodeData,

    # methods/monte_carlo/generic_ls_regression.jl
    generic_longstaff_schwartz_regression!,

    # quotes/Quotes.jl
    Quote,

    # currencies/currencies.jl
    NullCurrency, Currency,

    # InterestRates.jl
    ContinuousCompounding, SimpleCompounding, CompoundedCompounding, SimpleThenCompounded, ModifiedDuration,
    InterestRate, discount_factor, compound_factor, equivalent_rate, implied_rate,

    # exercise.jl
    AmericanExercise, BermudanExercise, EuropeanExercise,

    # Indexes
    IborIndex, LiborIndex, fixing_date, maturity_date, fixing, forecast_fixing, euribor_index, usd_libor_index, update_termstructure, is_valid_fixing_date, add_fixing!,

    # cash_flows/cash_flows.jl
    CouponMixin, accrual_start_date, accrual_end_date, ref_period_start, ref_period_end, SimpleCashFlow, Leg, ZeroCouponLeg, IRRFinder, operator, amount, date, duration, yield, previous_cashflow_date,
    accrual_days, accrual_days, next_cashflow, has_occurred, accrued_amount, next_coupon_rate, maturity_date, initialize!,

    # cash_flows/fixed_rate_coupon.jl
    FixedRateCoupon, FixedRateLeg,

    # cash_flows/floating_rate_coupon.jl
    BlackIborCouponPricer, IborCoupon, IborLeg, update_pricer!,

    # instruments/Instruments.jl
    update_pricing_engine, LongPosition, ShortPosition,

    # instruments/bond.jl
    BondMixin, FixedRateBond, FloatingRateBond, ZeroCouponBond, value, get_settlement_days, get_settlement_date, notional, accrued_amount, yield, duration, npv, clean_price, dirty_price, accrued_amount,

    # instruments/payoff.jl
    Put, Call, PlainVanillaPayoff,

    # instruments/option.jl
    VanillaOption,

    # instruments/forward.jl
    ForwardRateAgreement, forward_rate, spot_value, forward_value, implied_yield,

    # instruments/callability_schedule.jl
    DirtyCall, CleanCall, Price, Callability, CallabilitySchedule,

    # instruments/callable_bond.jl
    CallableFixedRateBond,

    # instruments/claim.jl
    FaceValueClaim,

    # instruments/swap.jl
    Payer, Receiver, SwapResults, VanillaSwap, NonstandardSwap, CreditDefaultSwap, fair_rate, update_all_ts!,

    # Swap Index
    SwapIndex, EuriborSwapIsdaFixA,

    # instruments/swaption.jl
    SettlementCash, SettlementPhysical, Swaption, NonstandardSwaption, calibration_basket,

    # termstructures/bond_helpers.jl
    FixedRateBondHelper, implied_quote,

    # termstructures/rate_helpers.jl
    SwapRateHelper, DepositRateHelper, FraRateHelper, implied_quote,

    # termstructures/TermStructures.jl
    check_range, max_date, time_from_reference,

    # termstructures/yield_term_structure.jl
    NullYieldTermStructure, FlatForwardTermStructure, JumpDate, JumpTime, InterpolatedDiscountCurve,
    calculated!, discount, zero_rate, forward_rate, discount_impl,

    # termstructures/bootstrap_traits.jl
    Discount, HazardRate, guess, min_value_after, max_value_after,

    # termstructures/curve.jl
    PiecewiseYieldCurve, PiecewiseDefaultCurve, FittedBondDiscountCurve, FittingCost, NullCurve,
    max_date, discount, calculate!, initialize!, value, nodes, survival_probability, default_probability, default_density, hazard_rate,

    # termstructures/vol_term_structure.jl
    ConstantOptionVolatility, ConstantSwaptionVolatility,

    # termstructures/black_vol_term_structure.jl
    BlackConstantVol,

    # termstructures/bootstrap.jl
    IterativeBootstrap, initialize, quote_error,

    # termstructures/credit_helper.jl
    SpreadCDSHelper,

    # termstructures/credit
    InterpolatedHazardRateCurve,

    # termstructures/nonlinear_fitting_methods.jl
    ExponentialSplinesFitting, SimplePolynomialFitting, NelsonSiegelFitting, SvenssonFitting, CubicBSplinesFitting, discount_function, guess_size,

    # Process
    # process/stochastic_process/jl
    OrnsteinUhlenbeckProcess, GsrProcess, expectation, variance,

    # process/black_scholes_process.jl
    BlackScholesMertonProcess,

    # process/heston_process.jl
    HestonProcess,

    #process/bates_process.jl
    BatesProcess,

    # models/parameter.jl
    ConstantParameter, G2FittingParameter, HullWhiteFittingParameter,

    # models/calibration_helpers.jl
    NaiveBasketType, SwaptionHelper, implied_volatility!, add_times_to!, model_value!, update_pricing_engine!, underlying_swap!, # for swaptionHelper only

    # models/equity/heston_model.jl
    HestonModel,

    # models/equity/bates_model.jl
    BatesModel,

    # models/short_rate/short_rate.jl
    PrivateConstraint, test, calibrate!, func_values, value, get_params,

    # models/short_rate/two_factor.jl
    G2,

    # models/short_rate/one_factor.jl
    BlackKarasinski, HullWhite, GSR, calibrate_volatilities_iterative!, get_volatilities,

    # models/market_models/callability/exercise_value.jl
    NothingExerciseValue,

    # models/market_models/brownian_generators.jl
    SobolDiagonalOrdering, SobolBrownianGeneratorFactory, create,

    # models/market_models/evolution_description.jl
    EvolutionDescription, money_market_measure,

    # models/market_models/products/multistep/multiproduct_multistep.jl
    MultiStepInverseFloater, ExerciseAdapter,

    # models/market_models/products/multistep/call_specified_multiproduct.jl
    CallSpecifiedMultiProduct,

    # models/market_models/products/pathwise/pathwise_product_calls_specified.jl
    CallSpecifiedPathwiseMultiProduct,

    # models/market_models/products/composite_product.jl
    MarketModelComposite, add_product!, finalize!,

    # models/market_models/correlations/exponential_correlation.jl
    ExponentialForwardCorrelation,

    # models/market_models/modesl/flatvol.jl
    FlatVol,

    # models/market_models/callability/swap_rate_trigger.jl
    SwapRateTrigger,

    # models/market_models/callability/market_model_basis_system.jl
    SwapForwardBasisSystem,

    # models/market_models/callability/collect_node_data.jl
    collect_node_data!,

    # models/market_models/callability/ls_strategy.jl
    LongstaffSchwartzExerciseStrategy,

    # models/market_models/pathwise_greeks/vega_bump_cluster.jl
    VegaBumpCollection,

    # models/market_models/pathwise_greeks/bump_instrument_jacobian.jl
    VolatilityBumpInstrumentJacobianSwaption, VolatilityBumpInstrumentJacobianCap, OrthogonalizedBumpFinder, get_vega_bumps!,

    # models/market_models/pathwise/pathwise_product_inverse_floater.jl
    MarketModelPathwiseInverseFloater,

    # models/market_models/evolvers/lognormal_fwd_rate_pc.jl
    LogNormalFwdRatePc,

    # models/market_models/evolvers/lognormal_fwd_rate_euler.jl
    LogNormalFwdRateEuler,

    # models/market_models/accounting_engine.jl
    AccountingEngine, multiple_path_values!,

    # models/market_models/pathwise_accounting_engine.jl
    PathwiseVegasOuterAccountingEngine, multiple_path_values!,

    # methods - finite difference
    FdmG2Solver, FdmHullWhiteSolver, CrankNelson, FdmSolverDesc, Fdm2DimSolver, FdmG2Op, FdmMesherComposite, FdmSimpleProcess1dMesher, FdmAffineModelSwapInnerValue,

    # methods - monte carlo
    MonteCarloModel, PathGenerator, add_samples!, Path,

    # pricing_engines/pricing_engines.jl
    DiscountingBondEngine, DiscountingSwapEngine, MidPointCdsEngine, TreeCallableFixedRateEngine,

    # pricing_engines/black_calculator.jl
    BlackCalculator, value, delta, vega,

    # pricing_engines/discretized_asset.jl
    DiscretizedSwaption, DiscretizedSwap,

    # pricing_engines/swaption_engines
    BlackSwaptionEngine, G2SwaptionEngine, JamshidianSwaptionEngine, TreeSwaptionEngine, FdG2SwaptionEngine, FdHullWhiteSwaptionEngine, Gaussian1DSwaptionEngine, Gaussian1DNonstandardSwaptionEngine,

    # pricing_engines/vanilla
    AnalyticEuropeanEngine, AnalyticHestonEngine, BatesEngine, BaroneAdesiWhaleyApproximationEngine, BjerksundStenslandApproximationEngine, IntegralEngine,
    FDEuropeanEngine, FDBermudanEngine, FDAmericanEngine, BinomialVanillaEngine, MCEuropeanEngine, MCAmericanEngine

# abstract types
include("abstract_types.jl")

# observer
include("observer.jl")

# lazy
include("lazy.jl")

# methods
include("methods/lattice.jl")
include("methods/bsm_lattice.jl")
include("methods/binomial_tree.jl")
include("methods/monte_carlo/utilities.jl")
include("methods/monte_carlo/brownian_bridge.jl")
include("methods/monte_carlo/generic_ls_regression.jl")

# Quotes ----------------------------
include("quotes/Quotes.jl")

# Currencies -----------------------
include("currencies/currencies.jl")

# Interest Rates ---------------------------------
include("InterestRates.jl")

# Exercise---------------------------------
include("exercise.jl")

# Indexes
include("indexes/indexes.jl")

# Cash Flows ------------------------------------
include("cash_flows/cash_flows.jl")
include("cash_flows/fixed_rate_coupon.jl")
include("cash_flows/floating_rate_coupon.jl")

# Instruments ------------------------
include("instruments/Instruments.jl")
# bond
include("instruments/bond.jl")
include("instruments/payoff.jl")
include("instruments/option.jl")
include("instruments/forward.jl")
include("instruments/callability_schedule.jl")
include("instruments/callable_bond.jl")
include("instruments/claim.jl")
include("instruments/swap.jl")
include("indexes/swap_index.jl") # need this here because swap index needs to know what a VanillaSwap is
include("instruments/swaption.jl")

# Term Structures -----------------------------------
include("termstructures/TermStructures.jl")
include("termstructures/curve.jl")
# helpers
include("termstructures/yield/bond_helpers.jl")
include("termstructures/yield/rate_helpers.jl")
include("termstructures/credit/credit_helpers.jl")
# bootstrapping
include("termstructures/bootstrap/bootstrap_traits.jl")
include("termstructures/bootstrap/bootstrap.jl")
# yield term structures
include("termstructures/yield/yield_term_structure.jl")
include("termstructures/yield/piecewise_yield_curve.jl")
include("termstructures/yield/fitted_bond_curve.jl")
include("termstructures/yield/nonlinear_fitting_methods.jl")
include("termstructures/yield/discount_curve.jl")
# volatility
include("termstructures/volatility/vol_term_structure.jl")
include("termstructures/volatility/black_vol_term_structure.jl")
include("termstructures/volatility/callable_bond_vol_structure.jl")
# credit
include("termstructures/credit/piecewise_default_curve.jl")
include("termstructures/credit/interpolated_hazard_rate_curve.jl")

# process
include("process/discretization.jl")
include("process/stochastic_process.jl")
include("process/black_scholes_process.jl")
include("process/heston_process.jl")
include("process/bates_process.jl")

# Models ---------------------------------
include("models/parameter.jl")
include("models/calibration_helpers.jl")
include("models/models.jl")
include("models/equity/heston_model.jl")
include("models/equity/bates_model.jl")
include("models/short_rate/short_rate.jl")
include("models/short_rate/two_factor.jl")
include("models/short_rate/one_factor/one_factor.jl")
include("models/short_rate/one_factor/black_karasinski.jl")
include("models/short_rate/one_factor/hull_white.jl")
include("models/short_rate/one_factor/gsr.jl")
include("models/market_models/utilities.jl")
include("models/market_models/evolution_description.jl")
include("models/market_models/brownian_generators.jl")
include("models/market_models/swap_forward_mappings.jl")
include("models/market_models/products/multistep/multiproduct_multistep.jl")
include("models/market_models/products/multistep/call_specified_multiproduct.jl")
include("models/market_models/products/composite_product.jl")
include("models/market_models/products/pathwise/pathwise_product_calls_specified.jl")
include("models/market_models/products/pathwise/pathwise_product_cashrebate.jl")
include("models/market_models/products/pathwise/rate_pseudo_root_jacobian.jl")
include("models/market_models/products/pathwise/pathwise_product_inverse_floater.jl")
include("models/market_models/models/flatvol.jl")
include("models/market_models/correlations/time_homogeneous_forward_correlation.jl")
include("models/market_models/correlations/exponential_correlation.jl")
include("models/market_models/callability/swap_rate_trigger.jl")
include("models/market_models/callability/exercise_value.jl")
include("models/market_models/callability/market_model_basis_system.jl")
include("models/market_models/callability/collect_node_data.jl")
include("models/market_models/pathwise_greeks/vega_bump_cluster.jl")
include("models/market_models/pathwise_greeks/swaption_pseudo_derivative.jl")
include("models/market_models/pathwise_greeks/bump_instrument_jacobian.jl")
include("models/market_models/callability/ls_strategy.jl")
include("models/market_models/curve_state.jl")
include("models/market_models/drift_computation.jl")
include("models/market_models/market_models.jl")
include("models/market_models/evolvers/lognormal_fwd_rate_euler.jl")
include("models/market_models/evolvers/lognormal_fwd_rate_pc.jl")
include("models/market_models/accounting_engine.jl")
include("models/market_models/pathwise_accounting_engine.jl")

# Finite Difference method
include("methods/finite_differences/fd_boundary_condition.jl")
include("methods/finite_differences/fd_pde.jl")
include("methods/finite_differences/fd_bsm_operator.jl")
include("methods/finite_differences/fd_layout.jl")
include("methods/finite_differences/fd_step_condition.jl")
include("methods/finite_differences/fd_mesher.jl")
include("methods/finite_differences/fd_calc.jl")
include("methods/finite_differences/fd_operator.jl")
include("methods/finite_differences/fd_scheme.jl")
include("methods/finite_differences/fd_model.jl")
include("methods/finite_differences/fd_solvers.jl")

# Monte Carlo method
include("methods/monte_carlo/path.jl")
include("methods/monte_carlo/sample.jl")
include("methods/monte_carlo/path_generator.jl")
include("methods/monte_carlo/lsm_basis_system.jl")
include("methods/monte_carlo/monte_carlo_model.jl")

# Pricing Engines ------------------------
include("pricing_engines/pricing_engines.jl")
include("pricing_engines/discretized_asset.jl")
include("pricing_engines/black_calculator.jl")
include("pricing_engines/black_formula.jl")
include("pricing_engines/mc_simulation.jl")
include("pricing_engines/lattice_shortrate_model_engine.jl")
include("pricing_engines/bond/discretized_bond.jl")
include("pricing_engines/bond/discounting_bond_engine.jl")
include("pricing_engines/bond/black_callable_bond_engine.jl")
include("pricing_engines/bond/tree_callable_bond_engine.jl")
include("pricing_engines/swap/discounting_swap_engine.jl")
include("pricing_engines/swap/discretized_swap.jl")
include("pricing_engines/credit/midpoint_cds_engine.jl")
include("pricing_engines/swaptions/swaption_engine.jl")
include("pricing_engines/swaptions/discretized_swaption.jl")
include("pricing_engines/swaptions/black_swaption_engine.jl")
include("pricing_engines/swaptions/G2_swaption_engine.jl")
include("pricing_engines/swaptions/fdg2_swaption_engine.jl")
include("pricing_engines/swaptions/fd_hullwhite_swaption_engine.jl")
include("pricing_engines/swaptions/jamshidian_swaption_engine.jl")
include("pricing_engines/swaptions/tree_swaption_engine.jl")
include("pricing_engines/swaptions/gaussian1d_swaption_engine.jl")
include("pricing_engines/swaptions/gaussian1d_nonstandard_swaption_engine.jl")
include("pricing_engines/vanilla/discretized_vanilla_option.jl")
include("pricing_engines/vanilla/analytic_european_engine.jl")
include("pricing_engines/vanilla/analytic_heston_engine.jl")
include("pricing_engines/vanilla/bates_engine.jl")
include("pricing_engines/vanilla/barone_adesi_whaley_engine.jl")
include("pricing_engines/vanilla/bjerksund_stensland_engine.jl")
include("pricing_engines/vanilla/integral_engine.jl")
include("pricing_engines/vanilla/fd_vanilla_engine.jl")
include("pricing_engines/vanilla/binomial_engine.jl")
include("pricing_engines/vanilla/mc_european_engine.jl")
include("pricing_engines/vanilla/mc_american_engine.jl")

# # Helpers NOW IN TERM STRUCTURE
# include("helpers/bond_helpers.jl")

type Settings
  evaluation_date::Date
end

settings = Settings(Date())

function set_eval_date!(sett::Settings, d::Date)
  sett.evaluation_date = d
end

export Settings, settings, set_eval_date!

end
