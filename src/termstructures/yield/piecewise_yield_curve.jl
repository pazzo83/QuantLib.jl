type PiecewiseYieldCurve{B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap} <: InterpolatedCurve{P, T}
  lazyMixin::LazyMixin
  settlementDays::Int
  referenceDate::Date
  instruments::Vector{B}
  dc::DC
  interp::P
  trait::T
  accuracy::Float64
  boot::BT
  times::Vector{Float64}
  dates::Vector{Date}
  data::Vector{Float64}
  errors::Vector{Function}
  validCurve::Bool
end

function PiecewiseYieldCurve{B <: BootstrapHelper, DC <: DayCount, P <: Interpolation, T <: BootstrapTrait, BT <: Bootstrap}(referenceDate::Date, instruments::Vector{B}, dc::DC, interp::P, trait::T,
                                accuracy::Float64, boot::BT = IterativeBootstrap())
  # get the initial length of instruments
  n = length(instruments)
  # create an initial state of the curve
  pyc = PiecewiseYieldCurve(LazyMixin(),
                            0,
                            referenceDate,
                            instruments,
                            dc,
                            interp,
                            trait,
                            accuracy,
                            boot,
                            Vector{Float64}(n + 1),
                            Vector{Date}(n + 1),
                            Vector{Float64}(n + 1),
                            Vector{Function}(n + 1),
                            false)

  # initialize the bootstrapping
  initialize(pyc.boot, pyc)

  return pyc
end
