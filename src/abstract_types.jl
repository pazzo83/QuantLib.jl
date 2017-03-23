# Observer
abstract type Observer end

# Lazy Object
abstract type LazyObject <: Observer end

# Process
abstract type StochasticProcess end
abstract type StochasticProcess1D <: StochasticProcess end
abstract type AbstractBlackScholesProcess <: StochasticProcess1D end
abstract type AbstractHestonProcess <: StochasticProcess1D end
abstract type AbstractDiscretization end
abstract type AbstractHestonDiscretization <: AbstractDiscretization end
abstract type BlackScholesType end

# Methods
abstract type Lattice end
abstract type TreeLattice <: Lattice end
abstract type AbstractTree end
abstract type AbstractBinomialTree <: AbstractTree end
abstract type BinomialTreeType end
abstract type EqualProbabilitiesBinomialTreeType <: BinomialTreeType end
abstract type EqualJumpsBinomialTreeType <: BinomialTreeType end
abstract type FdmSchemeDescType end
abstract type FdScheme end
abstract type FdmMesher end
abstract type Fdm1DMesher end
abstract type StepCondition end
abstract type CurveWrapper end
abstract type StepType end
abstract type FdmInnerValueCalculator end
abstract type FdmLinearOpComposite end
abstract type NinePointLinearOp end
abstract type BoundaryConditionType end
abstract type BC_Side end
# Monte Carlo
abstract type AbstractMonteCarloModel end
abstract type MCTrait end
abstract type AbstractPathPricer end
abstract type EarlyExercisePathPricer <: AbstractPathPricer end
abstract type LsmBasisSystemPolynomType end
abstract type LSMBasisSystemFunction <: Function end
# abstract BoundaryCondition

# Exercise
abstract type Exercise end
abstract type EarlyExercise <: Exercise end

# Instruments
abstract type Instrument <: LazyObject end
abstract type PositionType end
abstract type Bond <: Instrument end
abstract type AbstractCallableBond <: Bond end
abstract type AbstractRate <: Instrument end
abstract type Swap <: Instrument end
abstract type AbstractClaim end
abstract type SettlementType end
abstract type AbstractPayoff end
abstract type StrikedTypePayoff <: AbstractPayoff end
abstract type Option{E} <: Instrument end
abstract type OneAssetOption{E} <: Option{E} end
abstract type OptionType end
abstract type SwapType end
abstract type CDSProtectionSide end
abstract type Results end
abstract type CallType end
abstract type AbstractForward <: Instrument end

# Term Structures
abstract type TermStructure <: LazyObject end
abstract type YieldTermStructure <: TermStructure end
abstract type CreditTermStructure <: TermStructure end
# Curves
abstract type Curve <: YieldTermStructure end
abstract type InterpolatedCurve{P} <: Curve end
abstract type AbstractDefaultProbabilityTermStructure <: CreditTermStructure end
abstract type AbstractDefaultProbabilityCurve <: AbstractDefaultProbabilityTermStructure end
abstract type InterpolatedDefaultProbabilityCurve{P} <: AbstractDefaultProbabilityCurve end

abstract type VolatilityTermStructure <: TermStructure end
abstract type OptionletVolatilityStructure <: VolatilityTermStructure end
abstract type SwaptionVolatilityStructure <: VolatilityTermStructure end
abstract type BlackVolTermStructure <: VolatilityTermStructure end
abstract type LocalVolTermStructure <: VolatilityTermStructure end
abstract type CallableBondVolatilityStructure <: TermStructure end
abstract type AbstractSmileSection end
abstract type VolatilityType end
abstract type BootstrapTrait end
abstract type Bootstrap end
abstract type FittingMethod end
abstract type BootstrapHelper <: LazyObject end
abstract type BondHelper <: BootstrapHelper end
abstract type RateHelper <: BootstrapHelper end
abstract type AbstractCDSHelper <: BootstrapHelper end

# Pricing Engines
abstract type PricingEngine end
abstract type DiscretizedAsset end
abstract type DiscretizedOption <: DiscretizedAsset end
abstract type LatticeShortRateModelEngine{S} <: PricingEngine end
abstract type AbstractHestonEngine{HI} <: PricingEngine end
abstract type AbstractFDVanillaEngine <: PricingEngine end
abstract type FDMultiPeriodEngine <: AbstractFDVanillaEngine end
abstract type FDStepConditionEngine <: AbstractFDVanillaEngine end
abstract type AbstractVanillaEngine <: PricingEngine end
abstract type MCEngine{S, RSG} <: PricingEngine end
abstract type MCVanillaEngine{S, RSG} <: MCEngine{S, RSG} end
abstract type MCLongstaffSchwartzEngine{S, P, RSG} <: MCEngine{S, RSG} end
abstract type GaussianProbabilities end

# Cash Flows
abstract type CashFlows end
abstract type Leg <: CashFlows end
abstract type CashFlow end
abstract type Coupon <: CashFlow end
abstract type Duration end
abstract type CouponPricer end
abstract type IborCouponPricer <: CouponPricer end

# Indexes
abstract type InterestRateIndex end

# Models
abstract type Parameter end
abstract type CalibrationErrorType end
abstract type CalibrationHelper <: LazyObject end
abstract type CalibrationBasketType end
abstract type ModelType end
abstract type Model{T <: ModelType} <: LazyObject end
abstract type ShortRateModel{T} <: Model{T} end
abstract type OneFactorModel{T} <: ShortRateModel{T} end
abstract type Gaussian1DModel{T} <: OneFactorModel{T} end
abstract type TwoFactorModel{T} <: ShortRateModel{T} end
abstract type ComplexLogFormula end
abstract type HestonIntegration end
abstract type ShortRateDynamics end
abstract type ShortRateTree <: TreeLattice end
abstract type AbstractMarketModel end
abstract type AbstractMarketModelEvolver end
abstract type ExerciseStrategy end
abstract type MarketModelMultiProduct end
abstract type MarketModelBasisSystem end
abstract type MarketModelPathwiseMultiProduct end
abstract type MultiProductMultiStep <: MarketModelMultiProduct end
abstract type MarketModelExerciseValue end
abstract type PiecewiseConstantCorrelation end
abstract type BrownianGeneratorFactory end
abstract type BrownianGenerator end
abstract type SobolOrdering end
abstract type CurveState end

# Currencies
abstract type AbstractCurrency end

# Interest Rates
abstract type CompoundingType end
