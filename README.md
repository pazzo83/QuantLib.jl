# QuantLib.jl

[![Build Status](https://travis-ci.org/pazzo83/QuantLib.jl.svg?branch=master)](https://travis-ci.org/pazzo83/QuantLib.jl)

This package aims to provide a pure Julia version of the popular open source library QuantLib (written in C++ and interfaced with other languages via SWIG).  Right now the package is in an alpha state, but there is quite a bit of functionality already.

### Install
```julia
Pkg.clone("https://github.com/pazzo83/QuantLib.jl.git")
```

The package essentially contains the main QuantLib module and two sub-modules for various time-based and math-based operations.  Below is a fairly up-to-date status of what is included.

**Documentation**: <http://quantlibjl.readthedocs.org/en/latest/>

### Math
Interpolations:
* Backward Flat
* Linear
* Log Linear
* Cubic Spline
* BiCubic Spline (implemented with Dierckx)

Optimization methods:
* Simplex
* Levenberg Marquardt

Solvers:
* Brent
* Finite Differences
* Newton

### Time
Calendars (adopted from BusinessDays.jl and Ito.jl):
* Target (basically a null calendar with basic holidays)
* US Settlement Calendar
* US NYSE Calendar
* US NERC Calendar
* UK Settlement Calendar
* UK LSE Calendar
* UK LME Calendar

Day Counters:
* Actual 360
* Actual 365
* Bond Thirty 360
* Euro Bond Thirty 360
* ISMA Actual
* ISDA Actual
* AFBA Actual

### Instruments
Bonds:
* Fixed Rate Bond
* Floating Rate Bond

Options:
* Vanilla Option
* Swaption
* Nonstandard Swaption (used for Gaussian methods)

Swaps:
* Vanilla Swap
* Nonstandard Swap (used for Gaussian methods)
* Credit Default Swap (partial)

### Indexes
* Ibor
* Libor
* Euribor
* USD Libor
* Euribor Swap ISDA

### Methods
* Finite Differences
* Trinomial Tree
* Tree Lattice 1D & 2D
* Monte Carlo

### Models
Short Rate:
* Black Karasinski
* Gaussian Short Rate (GSR)
* Hull White
* G2

Equity:
* Bates Model
* Heston Model

Market Models
* Flat Vol

### Pricing Engines
Bond:
* Discounting Bond Engine
* Tree Callable Fixed Rate Bond Engine
* Black Callable Fixed Rate Bond Engine

Swap:
* Discounting Swap Engine

Credit:
* MidPoint CDS Engine

Swaptions:
* Black Swaption Engine
* Finite Differences Hull White Pricing Engine
* Finite Differences G2 Pricing Engine
* G2 Swaption Engine
* Gaussian 1D Nonstandard Swaption Engine
* Gaussian 1D Swaption Engine
* Jamshidian Swaption Engine
* Tree Swaption Engine

Vanilla:
* Analytic European Engine (for black scholes)
* Analytic Heston Engine
* Barone Adesi Whaley Engine
* Bates Engine
* Binomial Engine
* Bjerksund Stensland Approximation Engine
* FD Vanilla Engine
* Integral Engine
* MonteCarlo American Engine
* MonteCarlo European Engine

General:
* Black Scholes Calculator
* Black Formula
* MonteCarlo Simulation
* Lattice ShortRate Model Engine

### Processes
* Black Scholes Process
* Ornstein Uhlenbeck Process
* Gaussian Short Rate Process
* Bates Process
* Heston Process

### Term Structures
Credit:
* Piecewise Default Curve
* Interpolated Hazard Rate Curve

Volatility:
* Black Constant Vol
* Constant Optionlet Volatility
* Constant Swaption Volatility
* Local Constant Vol

Yield:
* Flat Forward
* Fitted Bond Curve (various fitting methods)
* Piecewise Yield Curve
* Discount Curve

## Example
### Price a fixed rate Bond
```julia
using QuantLib

settlement_date = Date(2008, 9, 18) # construct settlement date
# settings is a global singleton that contains global settings
set_eval_date!(settings, settlement_date - Dates.Day(3))

# settings that we will need to construct the yield curve
freq = QuantLib.Time.Semiannual()
tenor = QuantLib.Time.TenorPeriod(freq)
conv = QuantLib.Time.Unadjusted()
conv_depo = QuantLib.Time.ModifiedFollowing()
rule = QuantLib.Time.DateGenerationBackwards()
calendar = QuantLib.Time.USGovernmentBondCalendar()
dc_depo = QuantLib.Time.Actual365()
dc = QuantLib.Time.ISDAActualActual()
dc_bond = QuantLib.Time.ISMAActualActual()
fixing_days = 3

# build depos
depo_rates = [0.0096, 0.0145, 0.0194]
depo_tens = [Base.Dates.Month(3), Base.Dates.Month(6), Base.Dates.Month(12)]

# build bonds
issue_dates = [Date(2005, 3, 15), Date(2005, 6, 15), Date(2006, 6, 30), Date(2002, 11, 15),
              Date(1987, 5, 15)]
mat_dates = [Date(2010, 8, 31), Date(2011, 8, 31), Date(2013, 8, 31), Date(2018, 8, 15),
            Date(2038, 5, 15)]

coupon_rates = [0.02375, 0.04625, 0.03125, 0.04000, 0.04500]
market_quotes = [100.390625, 106.21875, 100.59375, 101.6875, 102.140625]

# construct the deposit and fixed rate bond helpers
insts = Vector{BootstrapHelper}(length(depo_rates) + length(issue_dates))
for i = 1:length(depo_rates)
  depo_quote = Quote(depo_rates[i])
  depo_tenor = QuantLib.Time.TenorPeriod(depo_tens[i])
  depo = DepositRateHelper(depo_quote, depo_tenor, fixing_days, calendar, conv_depo, true, dc_depo)
  insts[i] = depo
end

for i =1:length(coupon_rates)
  term_date = mat_dates[i]
  rate = coupon_rates[i]
  issue_date = issue_dates[i]
  market_quote = market_quotes[i]
  sched = QuantLib.Time.Schedule(issue_date, term_date, tenor, conv, conv, rule, true)
  bond = FixedRateBondHelper(Quote(market_quote), FixedRateBond(3, 100.0, sched, rate, dc_bond, conv,
                            100.0, issue_date, calendar, DiscountingBondEngine()))
  insts[i + length(depo_rates)] = bond
end

# Construct the Yield Curve
interp = QuantLib.Math.LogInterpolation()
trait = Discount()
bootstrap = IterativeBootstrap()
yts = PiecewiseYieldCurve(settlement_date, insts, dc, interp, trait, 0.00000000001, bootstrap)

# Build it
calculate!(yts)

# Build our Fixed Rate Bond
settlement_days = 3
face_amount = 100.0

fixed_schedule = QuantLib.Time.Schedule(Date(2007, 5, 15), Date(2017, 5, 15),
                QuantLib.Time.TenorPeriod(QuantLib.Time.Semiannual()), QuantLib.Time.Unadjusted(),
                QuantLib.Time.Unadjusted(), QuantLib.Time.DateGenerationBackwards(), false,
                QuantLib.Time.USGovernmentBondCalendar())

pe = DiscountingBondEngine(yts)

fixedrate_bond = FixedRateBond(settlement_days, face_amount, fixed_schedule, 0.045,
                  QuantLib.Time.ISMAActualActual(), QuantLib.Time.ModifiedFollowing(), 100.0,
                  Date(2007, 5, 15), fixed_schedule.cal, pe)

# Calculate NPV
npv(fixedrate_bond) # 107.66828913260542
```
