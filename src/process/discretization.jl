struct EulerDiscretization <: AbstractDiscretization end

drift(euler::EulerDiscretization, process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64) = drift(process, t0, x0) * dt

diffusion(euler::EulerDiscretization, process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64) = diffusion(process, t0, x0) * sqrt(dt)

function variance(euler::EulerDiscretization, process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64)
  sigma = diffusion(process, t0, x0)

  return sigma * sigma * dt
end
