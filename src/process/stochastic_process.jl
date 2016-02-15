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

lower_index(process::GsrProcess, x::Float64) = searchsortedlast(process.times, x) + 1
upper_index(process::GsrProcess, x::Float64) = searchsortedlast(process.times, x - 1e-15) + 2

function time2(process::GsrProcess, idx::Int)
  if idx == 1
    return 0.0
  end

  return process.times[idx - 1]
end

vol(process::GsrProcess, x::Int) = x > length(process.vols) ? process.vols[end] : process.vols[x]
rev_zero(process::GsrProcess, x::Int) = x > length(process.revZero) ? process.revZero[end] : process.revZero[x]
rev(process::GsrProcess, x::Int) = x > length(process.reversions) ? process.reversions[end] : process.reversions[x]

floored_time(process::GsrProcess, x::Int, flr::Float64) = floor != -1.0 ? max(floor, time2(process, index)) : time2(process, index)
capped_time(process::GsrProcess, x::Int, cap::Float64) = cap != -1.0 ? min(cap, time2(process, index)) : time2(process, index)

function variance!(process::GsrProcess, w::Float64, ::Float64, dt::Float64)
  t = w + dt
  key = Pair(w, t)

  p = get(cache3, key, typemax(Float64))
  if p != typemax(Float64)
    return p
  end

  res = 0.0
  for k = lower_index(process, w):upper_index(process, t) - 1
    res2 = vol(process, k) * process(vol, k)
    res2 *= rev_zero(process, k) ? -(floored_time(process, k, w) - capped_time(process, k + 1, t)) :
                                  (1.0 - exp(2.0 * rev(process, k) * (floored_time(process, k, w) - capped_time(process, k + 1, t)))) /
                                  (2.0 * rev(process, k))

    for i = k + 1:upper_index(process, t) - 1
      res2 *= exp(-2.0 * rev(process, i) * (capped_time(process, i + 1, t) - time2(process, i)))
    end

    res += res2
  end

  process.cache3[key] = res

  return res
end

function expectationp1(process::GsrProcess, w::Float64, xw::Float64, dt::Float64)
  t = w + dt
  key = Pair(w, t)
  k = get(process.cache1, key, typemax(Float64))
  if k != typemax(Float64)
    return xw * k
  end

  res2 = 1.0
  for i = lower_index(process, w):upper_index(process, t) - 1
    res2 *= exp(-rev(process, i) * (capped_time(process, i + 1, t) - floored_time(process, i, w)))
  end

  process.cache1[key] = res2
  return res2 * xw
end

function expectationp2(process::GsrProcess, w::Float64, dt::Float64)
  t = w * dt
  key = Pair(w, t)

  p = get(process.cache2, key, typemax(Float64))
  if p != typemax(Float64)
    return p
  end

  T = process.T

  res = 0.0
  for k = lower_index(process, w):upper_index(process, t) - 1
    for l = 1:k - 1
      res2 = 1.0
      # alpha 1
      res2 *= rev_zero(process, l) ? vol(process, l) * vol(process, l) * (time2(process, l + 1) - time2(process, l)) :
                                   vol(process, l) * vol(process, l) / (2.0 * rev(process, l)) * (1.0 - exp(-2.0 * rev(process, l) *
                                                                          (time2(process, l + 1) - time2(process, l))))
      # zeta_i (i > k)
      for i = k+1:upper_index(i) - 1
        res2 *= exp(-rev(process, i) * (capped_time(process, i + 1, t) - time2(process, i)))
      end
      # beta_i (j<k)
      for j = l+1:k - 1
        res2 *= exp(-2.0 * rev(process, j) * (time2(process, j + 1) - time2(process, j)))
      end

      # zeta_k beta_k
    end
    res2 = 1.0
    # alpha_k zeta_k
    res2 *= rev_zero(process, k) ? vol(process, k) * vol(process, k) / 4.0 * (4.0 * ^(capped_time(process, k+1, t) - time2(process, k), 2.0) -
                                    (^(floored_time(process, k, w) - 2.0 * time2(process, k) + capped_time(process, k+1, t), 2.0) +
                                    ^(capped_time(process, k+1, t) - floored_time(process, k, w), 2.0))) :
                                 vol(process, k) * vol(process, k) / (2.0 * rev(process, k) * rev(process, k)) * (exp(-2.0 * rev(process, k) *
                                    (capped_time(process, k+1, t) - time2(process, k))) + 1.0 - (exp(-rev(process, k) * (floored_time(process, k, w) -
                                    2.0 * time2(process, k) + capped_time(process, k+1, t))) + exp(-rev(process, k) * (capped_time(process, k+1, t) -
                                    floored_time(process, k, w)))))

    # zeta_i (i>k)
    for i = k+1:upper_index(process, t) - 1
      res2 *= exp(-rev(process, i) * (capped_time(process, i+1, t) - time2(process, i)))
    end
    res += res2
  end
end

function expectation(process::GsrProcess, w::Float64, xw::Float64, dt::Float64)
  # t = w + dt

  return expectationp1(process, w, xw, dt) + expectationp2(process, w, dt)
end

std_deviation(process::GsrProcess, t0::Float64, x0::Float64, dt::Float64) = sqrt(variance(process, t0, x0, dt))

evolve(process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64, dw::Float64) =
      apply(process, expectation(process, t0, x0, dt), std_deviation(process, t0, x0, dt) * dw)

apply(process::StochasticProcess1D, x0::Float64, dx::Float64) = x0 + dx
