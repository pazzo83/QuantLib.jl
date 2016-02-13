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

type GsrProcess <: StochasticProcess1D
  times::Vector{Float64}
  vols::Vector{Float64}
  reversions::Vector{Float64}
  T::Float64
  revZero::Vector{Bool}
  cache1::Dict{Pair{Float64, Float64}, Float64}
  cache2::Dict{Pair{Float64, Float64}, Float64}
  cache3::Dict{Pair{Float64, Float64}, Float64}
  cache4::Dict{Float64, Float64}
  cache5::Dict{Pair{Float64, Float64}, Float64}
end

function GsrProcess(times::Vector{Float64}, vols::Vector{Float64}, reversions::Vector{Float64}, T::Float64 = 60.0)
  revZero = falses(length(reversions))
  cache1 = Dict{Pair{Float64, Float64}, Float64}()
  cache2 = Dict{Pair{Float64, Float64}, Float64}()
  cache3 = Dict{Pair{Float64, Float64}, Float64}()
  cache4 = Dict{Float64, Float64}()
  cache5 = Dict{Pair{Float64, Float64}, Float64}()

  return GsrProcess(times, vols, reversions, T, revZero, cache1, cache2, cache3, cache4, cache5)
end

function flush_cache!(gsrP::GsrProcess)
  for i in eachindex(gsrP.reversions)
    if abs(gsrP.reversions[i]) < 1e-4
      gsrP.revZero[i] = true
    else
      gsrP.revZero[i] = false
    end
  end
  empty!(gsrP.cache1)
  empty!(gsrP.cache2)
  empty!(gsrP.cache3)
  empty!(gsrP.cache4)
  empty!(gsrP.cache4)

  return gsrP
end

evolve(process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64, dw::Float64) =
      apply(process, expectation(process, t0, x0, dt), std_deviation(process, t0, x0, dt) * dw)

apply(process::StochasticProcess1D, x0::Float64, dx::Float64) = x0 + dx
