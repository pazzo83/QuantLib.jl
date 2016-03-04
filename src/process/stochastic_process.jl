type OrnsteinUhlenbeckProcess <: StochasticProcess1D
  speed::Float64
  vol::Float64
  x0::Float64
  level::Float64

  OrnsteinUhlenbeckProcess(speed::Float64, vol::Float64, x0::Float64 = 0.0, level::Float64 = 0.0) = new(speed, vol, x0, level)
end

expectation(process::OrnsteinUhlenbeckProcess, ::Float64, x0::Float64, dt::Float64) = process.level + (x0 - process.level) * exp(-process.speed * dt)
get_x0(process::OrnsteinUhlenbeckProcess) = process.x0

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

get_x0(::GsrProcess) = 0.0

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

lower_index(process::GsrProcess, x::Float64) = upper_bound(process.times, x)
upper_index(process::GsrProcess, x::Float64) = upper_bound(process.times, x - 1e-15) + 1

function time2(process::GsrProcess, idx::Int)
  if idx == 1
    return 0.0
  end

  if idx > length(process.times) + 1
    return process.T
  end

  return process.times[idx - 1]
end

vol(process::GsrProcess, x::Int) = x > length(process.vols) ? process.vols[end] : process.vols[x]
rev_zero(process::GsrProcess, x::Int) = x > length(process.revZero) ? process.revZero[end] : process.revZero[x]
rev(process::GsrProcess, x::Int) = x > length(process.reversions) ? process.reversions[end] : process.reversions[x]

floored_time(process::GsrProcess, x::Int, flr::Float64) = flr != -1.0 ? max(flr, time2(process, x)) : time2(process, x)
capped_time(process::GsrProcess, x::Int, cap::Float64) = cap != -1.0 ? min(cap, time2(process, x)) : time2(process, x)

function variance(process::GsrProcess, w::Float64, ::Float64, dt::Float64)
  t = w + dt
  key = Pair(w, t)

  p = get(process.cache3, key, typemax(Float64))
  if p != typemax(Float64)
    return p
  end

  res = 0.0
  for k = lower_index(process, w):upper_index(process, t) - 1
    res2 = vol(process, k) * vol(process, k)
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
  t = w + dt
  key = Pair(w, t)

  p = get(process.cache2, key, typemax(Float64))
  if p != typemax(Float64)
    return p
  end

  T = process.T

  res = 0.0
  # int A(s, t)y(s)
  for k = lower_index(process, w):upper_index(process, t) - 1
    # l < k
    for l = 1:k - 1
      res2 = 1.0
      # alpha 1
      if rev_zero(process, l)
        res2 *= vol(process, l) * vol(process, l) * (time2(process, l + 1) - time2(process, l))
      else
        res2 *= vol(process, l) * vol(process, l) / (2.0 * rev(process, l)) * (1.0 - exp(-2.0 * rev(process, l) *
                (time2(process, l + 1) - time2(process, l))))
      end
      # zeta_i (i > k)
      for i = k+1:upper_index(process, t) - 1
        res2 *= exp(-rev(process, i) * (capped_time(process, i + 1, t) - time2(process, i)))
      end
      # beta_i (j<k)
      for j = l+1:k - 1
        res2 *= exp(-2.0 * rev(process, j) * (time2(process, j + 1) - time2(process, j)))
      end

      # zeta_k beta_k
      if rev_zero(process, k)
        res2 *= 2.0 * time2(process, k) - floored_time(process, k, w) - capped_time(process, k + 1, t) -
                2.0 * (time2(process, k) - capped_time(process, k + 1, t))
      else
        res2 *= (exp(rev(process, k) * (2.0 * time2(process, k) - floored_time(process, k, w) - capped_time(process, k + 1, t))) -
                exp(2.0 * rev(process, k) * (time2(process, k) - capped_time(process, k + 1, t)))) / rev(process, k)
      end
      res += res2
    end
    # l = k
    res2 = 1.0
    # alpha_k zeta_k
    if rev_zero(process, k)
      res2 *= vol(process, k) * vol(process, k) / 4.0 * (4.0 * ^(capped_time(process, k + 1, t) - time2(process, k), 2.0) -
              (^(floored_time(process, k, w) - 2.0 * time2(process, k) + capped_time(process, k + 1, t), 2.0) +
              ^(capped_time(process, k + 1, t) - floored_time(process, k, w), 2.0)))
    else
      res2 *= vol(process, k) * vol(process, k) / (2.0 * rev(process, k) * rev(process, k)) * (exp(-2.0 * rev(process, k) *
              (capped_time(process, k + 1, t) - time2(process, k))) + 1.0 - (exp(-rev(process, k) * (floored_time(process, k, w) -
              2.0 * time2(process, k) + capped_time(process, k + 1, t))) + exp(-rev(process, k) *
              (capped_time(process, k + 1, t) - floored_time(process, k, w)))))
    end

    # zeta_i (i>k)
    for i = k+1:upper_index(process, t) - 1
      res2 *= exp(-rev(process, i) * (capped_time(process, i + 1, t) - time2(process, i)))
    end
    # no beta_j in this case
    res += res2
  end

  # int -A(s, t) sigma^2 G(s, T)
  for k = lower_index(process, w):upper_index(process, t) - 1
    res2 = 0.0
    # l > k
    for l = k+1:upper_index(process, T) - 1
      res3 = 1.0
      # eta_l
      res3 *= rev_zero(process, l) ? capped_time(process, l + 1, T) - time2(process, l) :
                                      (1.0 - exp(-rev(process, l) * (capped_time(process, l + 1, T) - time2(process, l)))) /
                                      rev(process, l)

      # zeta_i (i > k)
      for i = k+1:upper_index(process, t) - 1
        res3 *= exp(-rev(process, i) * (capped_time(process, i + 1, t) - time2(process, i)))
      end

      # gamma_j (j > k)
      for j = k+1:l-1
        res3 *= exp(-rev(process, j) * (time2(process, j + 1) - time2(process, j)))
      end

      # zeta_k gamma_k
      if rev_zero(process, k)
        res3 *= (capped_time(process, k + 1, t) - time2(process, k + 1) - (2.0 * floored_time(process, k, w) - capped_time(process, k + 1, t) -
                time2(process, k + 1))) / 2.0
      else
        res3 *= (exp(rev(process, k) * (capped_time(process, k + 1, t) - time2(process, k + 1))) - exp(rev(process, k) *
                (2.0 * floored_time(process, k, w) - capped_time(process, k + 1, t) - time2(process, k + 1)))) / (2.0 * rev(process, k))
      end

      res2 += res3
    end
    # l = k
    res3 = 1.0
    # eta_k zeta_k
    if rev_zero(process, k)
      res3 *= (-^(capped_time(process, k + 1, t) - capped_time(process, k + 1, T), 2.0) - 2.0 * ^(capped_time(process, k + 1, t) - floored_time(process, k, w), 2.0) +
              ^(2.0 * floored_time(process, k, w) - capped_time(process, k + 1, T) - capped_time(process, k + 1, t), 2.0)) / 4.0
    else
      res3 *= (2.0 - exp(rev(process, k) * (capped_time(process, k + 1, t) - capped_time(process, k + 1, T))) -
              (2.0 * exp(-rev(process, k) * (capped_time(process, k + 1, t) - floored_time(process, k, w))) -
              exp(rev(process, k) * (2.0 * floored_time(process, k, w) - capped_time(process, k + 1, T) -
              capped_time(process, k + 1, t))))) / (2.0 * rev(process, k) * rev(process, k))
    end

    # zeta_i (i > k)
    for i = k+1:upper_index(process, t) - 1
      res3 *= exp(-rev(process, i) * (capped_time(process, i + 1, t) - time2(process, i)))
    end

    # no gamma_j in this case
    res2 += res3

    # add to main accumulator
    res += -vol(process, k) * vol(process, k) * res2
  end
  process.cache2[key] = res

  return res
end

function G!(process::GsrProcess, t::Float64, w::Float64, ::Float64)
  key = Pair(w, t)
  k = get(process.cache5, key, typemax(Float64))
  if k != typemax(Float64)
    return k
  end

  res = 0.0
  for i = lower_index(process, t):upper_index(process, w) - 1
    res2 = 1.0
    for j = lower_index(process, t):i - 1
      res2 *= exp(-rev(process, j) * (time2(process, j + 1) - floored_time(process, j, t)))
    end
    res2 *= rev_zero(process, i) ? capped_time(process, i + 1, w) - floored_time(process, i, t) :
                                  (1.0 - exp(-rev(process, i) * (capped_time(process, i + 1, w) - floored_time(process, i, t)))) / rev(process, i)

    res += res2
  end

  process.cache5[key] = res
  return res
end

function y!(process::GsrProcess, t::Float64)
  key = t
  k = get(process.cache4, key, typemax(Float64))
  if k != typemax(Float64)
    return k
  end

  res = 0.0
  for i = 1:upper_index(process, t) - 1
    res2 = 1.0
    for j = i + 1:upper_index(process, t) - 1
      res2 *= exp(-2.0 * rev(process, j) * (capped_time(process, j + 1, t) - time2(process, j)))
    end
    res2 *= rev_zero(process, i) ? vol(process, i) * vol(process, i) * (capped_time(process, i+1, t) - time2(process, i)) :
                                   (vol(process, i) * vol(process, i) / (2.0 * rev(process, i)) * (1.0 - exp(-2.0 * rev(process, i) * (capped_time(process, i+1, t) - time2(process, i)))))

    res += res2
  end

  process.cache4[key] = res
  return res
end

function expectation(process::GsrProcess, w::Float64, xw::Float64, dt::Float64)
  # t = w + dt

  # p1 = expectationp1(process, w, xw, dt)
  # p2 = expectationp2(process, w, dt)
  # println(p1)
  # println(p2)
  # error("DIE")

  return expectationp1(process, w, xw, dt) + expectationp2(process, w, dt)
end

std_deviation(process::GsrProcess, t0::Float64, x0::Float64, dt::Float64) = sqrt(variance(process, t0, x0, dt))

evolve(process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64, dw::Float64) =
      apply(process, expectation(process, t0, x0, dt), std_deviation(process, t0, x0, dt) * dw)

apply(process::StochasticProcess1D, x0::Float64, dx::Float64) = x0 + dx

std_deviation(process::StochasticProcess1D, t0::Float64, x0::Float64, dt::Float64) = diffusion(process.disc, process, t0, x0, dt)
