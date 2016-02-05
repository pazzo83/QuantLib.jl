type OrnsteinUhlenbeckProcess <: StochasticProcess1D
  speed::Float64
  vol::Float64
  x0::Float64
  level::Float64

  OrnsteinUhlenbeckProcess(speed::Float64, vol::Float64, x0::Float64 = 0.0, level::Float64 = 0.0) = new(speed, vol, x0, level)
end

expectation(process::OrnsteinUhlenbeckProcess, ::Float64, x0::Float64, dt::Float64) = process.level + (x0 - process.level) * exp(-process.speed * dt)

function variance(process::OrnsteinUhlenbeckProcess, ::Float64, ::Float64, dt::Float64)
  v = process.vol
  if process.speed < sqrt(eps())
    # algebraic limits for small speed
    return v * v * dt
  else
    return 0.5 * v * v / process.speed * (1.0 - exp(-2.0 * process.speed * dt))
  end
end

std_deviation(process::OrnsteinUhlenbeckProcess, t::Float64, x0::Float64, dt::Float64) = sqrt(variance(process, t, x0, dt))

evolve(process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64, dw::Float64) =
      apply(process, expectation(process, t0, x0, dt), std_deviation(process, t0, x0, dt) * dw)

apply(process::StochasticProcess1D, x0::Float64, dx::Float64) = x0 + dx
