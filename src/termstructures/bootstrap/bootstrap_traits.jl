# bootstrapping traits
const avg_rate = 0.05
const max_rate = 1.0

const AVG_HAZARD_RATE = 0.01
const MAX_HAZARD_RATE = 1.0

# Discount bootstrap trait
type Discount <: BootstrapTrait end
type HazardRate <: BootstrapTrait end

# Type aliases
typealias InterpolatedHazardRateCurve{P, T} InterpolatedDefaultProbabilityCurve{P, HazardRate}

initial_value(::Discount) = 1.0
max_iterations(::Discount) = 100

initial_value(::HazardRate) = AVG_HAZARD_RATE
max_iterations(::HazardRate) = 30

function guess{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, valid::Bool)
  if valid
    # return previous iteration value
    return ts.data[i]
  end

  if i == 2
    # first pillar
    return 1.0 / (1.0 + avg_rate * ts.times[2])
  end

  r = -log(ts.data[i - 1]) / ts.times[i - 1]
  return exp(-r * ts.times[i])
end

function min_value_after{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, valid::Bool)
  if valid
    return ts.data[end] / 2.0
  end

  dt = ts.times[i] - ts.times[i - 1]
  return ts.data[i - 1] * exp( -max_rate * dt)
end

function max_value_after{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, ::Bool)
  dt = ts.times[i] - ts.times[i - 1]
  return ts.data[i - 1] * exp(max_rate * dt)
  #return ts.data[i - 1]
end

function update_guess!{I <: Integer, T <: TermStructure}(::Discount, i::I, ts::T, discount::Float64)
  ts.data[i] = discount
  return ts
end

function guess(::HazardRate, i::Int, ts::TermStructure, valid::Bool)
  if valid
    return ts.data[i]
  end

  if i == 2 # first pillar
    return AVG_HAZARD_RATE
  end

  # extrapolate
  d = ts.dates[i]
  return hazard_rate(ts, d)
end

function min_value_after(::HazardRate, i::Int, ts::TermStructure, valid::Bool)
  if valid
    r = min(ts.data[1], ts.data[end])
    return r / 2.0
  end

  return eps()
end

function max_value_after(::HazardRate, i::Int, ts::TermStructure, valid::Bool)
  if valid
    r = max(ts.data[1], ts.data[end])
    return r * 2.0
  end

  return MAX_HAZARD_RATE
end

function update_guess!(::HazardRate, i::Int, ts::TermStructure, rate::Float64)
  ts.data[i] = rate
  if i == 2
    ts.data[1] = rate # first point handled as well
  end

  return ts
end
