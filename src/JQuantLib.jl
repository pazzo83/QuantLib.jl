# JQuantLib
module JQuantLib

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

# Time module
include("time/Time.jl")

# Math module
include("math/Math.jl")

# MAIN MODULE CODE
using JQuantLib.Math, JQuantLib.Time

# Constants #
const basisPoint = 0.0001

export
    # abstract_types.jl
    LazyObject,

    Exercise, EarlyExercise, CompoundingType, TermStructure, YieldTermStructure, InterpolatedCurve, BootstrapTrait, Bootstrap, BootstrapHelper, BondHelper, RateHelper,
    FittingMethod, CashFlows, CashFlow, Coupon, CouponPricer, IborCouponPricer, Instrument, Bond, Swap, SwapType, PricingEngine, Duration, AbstractRate, Results,
    InterestRateIndex, AbstractCurrency, Parameter, CalibrationHelper, ShortRateModel, FdmMesher,

    # lazy.jl
    LazyMixin, calculate!, recalculate!,

    # Process
    # process/stochastic_process/jl
    OrnsteinUhlenbeckProcess, expectation, variance,

    # Methods
    # methods/lattice.jl
    TreeLattice1D, Branching, TrinomialTree,

    # quotes/Quotes.jl
    Quote,

    # currencies/currencies.jl
    NullCurrency, Currency,

    # InterestRates.jl
    ContinuousCompounding, SimpleCompounding, CompoundedCompounding, SimpleThenCompounded, ModifiedDuration,
    InterestRate, discount_factor, compound_factor, equivalent_rate, implied_rate,

    # exericise.jl
    AmericanExercise, BermudanExercise, EuropeanExercise,

    # Indexes
    IborIndex, LiborIndex, fixing_date, maturity_date, fixing, forecast_fixing, euribor_index, usd_libor_index,

    # cash_flows/cash_flows.jl
    CouponMixin, accrual_start_date, accrual_end_date, ref_period_start, ref_period_end, SimpleCashFlow, Leg, ZeroCouponLeg, IRRFinder, operator, amount, date, duration, yield, previous_cashflow_date,
    accrual_days, accrual_days, next_cashflow, has_occurred, next_coupon_rate, initialize!,

    # cash_flows/fixed_rate_coupon.jl
    FixedRateCoupon, FixedRateLeg,

    # cash_flows/floating_rate_coupon.jl
    BlackIborCouponPricer, IborCoupon, IborLeg, update_pricer!,

    # instruments/Instruments.jl
    update_pricing_engine,

    # instruments/bond.jl
    BondMixin, FixedRateBond, FloatingRateBond, ZeroCouponBond, value, get_settlement_days, get_settlement_date, notional, accrued_amount, yield, duration, npv, clean_price, dirty_price, accrued_amount,

    # instruments/option.jl
    Put, Call,

    # instruments/claim.jl
    FaceValueClaim,

    # instruments/swap.jl
    Payer, Receiver, SwapResults, VanillaSwap, CreditDefaultSwap, fair_rate,

    # instruments/swaption.jl
    SettlementCash, SettlementPhysical, Swaption,

    # termstructures/bond_helpers.jl
    FixedRateBondHelper, implied_quote,

    # termstructures/rate_helpers.jl
    SwapRateHelper, DepositRateHelper, implied_quote,

    # termstructures/TermStructures.jl
    check_range, max_date, time_from_reference,

    # termstructures/yield_term_structure.jl
    NullYieldTermStructure, FlatForwardTermStructure, JumpDate, JumpTime,
    calculated!, discount, zero_rate, forward_rate, discount_impl,

    # termstructures/bootstrap_traits.jl
    Discount, HazardRate, guess, min_value_after, max_value_after,

    # termstructures/curve.jl
    PiecewiseYieldCurve, PiecewiseDefaultCurve, FittedBondDiscountCurve, FittingCost, NullCurve,
    max_date, discount, calculate!, initialize!, value, nodes, survival_probability, 

    # termstructures/vol_term_structure.jl
    ConstantOptionVolatility, ConstantSwaptionVolatility,

    # termstructures/bootstrap.jl
    IterativeBootstrap, initialize, quote_error,

    # termstructures/credit_helper.jl
    SpreadCDSHelper,

    # termstructures/nonlinear_fitting_methods.jl
    ExponentialSplinesFitting, SimplePolynomialFitting, NelsonSiegelFitting, SvenssonFitting, CubicBSplinesFitting, discount_function, guess_size,

    # models/parameter.jl
    ConstantParameter, G2FittingParameter, HullWhiteFittingParameter,

    # models/calibration_helpers.jl
    SwaptionHelper, implied_volatility!, add_times_to!, model_value!, update_pricing_engine!, # for swaptionHelper only

    # models/short_rate/short_rate.jl
    PrivateConstraint, test, calibrate!, func_values, value, get_params,

    # models/short_rate/two_factor.jl
    G2,

    # models/short_rate/one_factor.jl
    BlackKarasinski, HullWhite,

    # methods - finite difference
    FdmG2Solver,FdmHullWhiteSolver,

    # pricing_engines/pricing_engines.jl
    DiscountingBondEngine, DiscountingSwapEngine, MidPointCdsEngine,

    # pricing_engines/discretized_asset.jl
    DiscretizedSwaption, DiscretizedSwap,

    # pricing_engines/swaption_engines.jl
    G2SwaptionEngine, JamshidianSwaptionEngine, TreeSwaptionEngine, FdG2SwaptionEngine, FdHullWhiteSwaptionEngine

# abstract types
include("abstract_types.jl")

# observer
include("observer.jl")

# lazy
include("lazy.jl")

# process
include("process/stochastic_process.jl")

# methods
include("methods/lattice.jl")

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
include("instruments/option.jl")
include("instruments/claim.jl")
include("instruments/swap.jl")
include("instruments/swaption.jl")

# helpers
include("termstructures/bond_helpers.jl")
include("termstructures/rate_helpers.jl")

# Term Structures -----------------------------------
include("termstructures/TermStructures.jl")
# yield term structures
include("termstructures/yield_term_structure.jl")
# bootstrapping traits
include("termstructures/bootstrap_traits.jl")
# Curves
include("termstructures/curve.jl")
# volatility
include("termstructures/vol_term_structure.jl")
# bootstrapping
include("termstructures/bootstrap.jl")
# credit helpers
include("termstructures/credit_helpers.jl")

# nonlinear fitting methods
include("termstructures/nonlinear_fitting_methods.jl")

# Models ---------------------------------
include("models/parameter.jl")
include("models/calibration_helpers.jl")
include("models/models.jl")
include("models/short_rate/short_rate.jl")
include("models/short_rate/two_factor.jl")
include("models/short_rate/one_factor.jl")

# Finite Difference method
include("methods/finite_differences/fd_layout.jl")
include("methods/finite_differences/fd_step_condition.jl")
include("methods/finite_differences/fd_mesher.jl")
include("methods/finite_differences/fd_calc.jl")
include("methods/finite_differences/fd_operator.jl")
include("methods/finite_differences/fd_scheme.jl")
include("methods/finite_differences/fd_model.jl")
include("methods/finite_differences/fd_solvers.jl")

# Pricing Engines ------------------------
include("pricing_engines/pricing_engines.jl")
include("pricing_engines/discretized_asset.jl")
include("pricing_engines/swaption_engines.jl")

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
