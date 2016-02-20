# JQuantLib

[![Build Status](https://travis-ci.org/pazzo83/JQuantLib.jl.svg?branch=master)](https://travis-ci.org/pazzo83/JQuantLib.jl)

This package aims to provide a pure Julia version of the popular open source library QuantLib (written in C++ and interfaced with other languages via SWIG).  Right now the package is in an alpha state, but there is quite a bit of functionality already.

The package is essentially contains the main JQuantLib module and two sub-modules for various time-based and math-based operations.  Below is a fairly up-to-date status of what is included.

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

### Models
Short Rate:
* Black Karasinski
* Gaussian Short Rate (GSR)
* Hull White
* G2

### Pricing Engines
Bond:
* Discounting Bond Engine

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

General:
* Black Scholes Calculator

### Processes
* Black Scholes Process
* Ornstein Uhlenbeck Process
* Gaussian Short Rate Process

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
