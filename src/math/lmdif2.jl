const MACHEP = 1.2e-16
const DWARF = 1.0e-38

function jacFcnDefault!(m::Int, n::Int, x::Vector{Float64}, fvec::Vector{Float64}, fjac::Matrix{Float64}, eps_::Float64, wa::Vector{Float64}, fcn!::Function)
  for j = 1:n
    temp = x[j]
    step = max(eps_ * eps_, eps_ * abs(temp))

    x[j] += step
    fcn!(m, n, x, wa)
    for i = 1:m
      fjac[i, j] = (wa[i] - fvec[i]) / step
    end

    #restore
    x[j] = temp
  end
  return fjac
end

function lm_qrsolv!(n::Int, r::Matrix{Float64}, ldr::Int, pivot::Vector{Int}, diagonal::Vector{Float64}, qtb::Vector{Float64}, x::Vector{Float64}, sdiag::Vector{Float64}, W::Vector{Float64})

  # Copy R and Q'* b to preserve input and initialize S
  # In particular, save the diagonal elements of R in x

  # instead of loop...
  # r = triu(r) + triu(r)'
  # x = diag(r)
  # W = qtb

  for j = 1:n
    for i = j:n
      r[(j-1) * ldr + i] = r[(i-1)*ldr + j]
    end

    x[j] = r[(j-1) * ldr + j]
    W[j] = qtb[j]
  end

  # eliminate the diagonal matrix D using a Givens rotation
  for j = 1:n
    if diagonal[pivot[j]] != 0.0
      for k=j:n
        sdiag[k] = 0.0
      end
      sdiag[j] = diagonal[pivot[j]]

      # the transformation to elminate the row of D modify only a single element of Q' * b beyond the first n, which is initially 0
      qtbpj = 0.0
      for k = j:n
        # determine the givens rotation which elminates the appropriate element in the current row of D
        if sdiag[k] == 0.0
          continue
        end

        kk = k + ldr * (k - 1)
        if abs(r[kk]) < abs(sdiag[k])
          _cot = r[kk] / sdiag[k]
          _sin = 1.0 / hypot(1.0, _cot)
          _cos = _sin * _cot
        else
          _tan = sdiag[k] / r[kk]
          _cos = 1.0 / hypot(1.0, _tan)
          _sin = _cos * _tan
        end

        # compute the modified diagonal element of R and the modified element of Q' * b, 0
        r[kk] = _cos * r[kk] + _sin * sdiag[k]
        temp = _cos * W[k] + _sin * qtbpj
        qtbpj = -_sin * W[k] + _cos * qtbpj
        W[k] = temp

        # accumulate the transformation in the row of S
        for i = k+1:n
          temp = _cos * r[(k - 1) * ldr + i] + _sin * sdiag[i]
          sdiag[i] = -_sin * r[(k-1) * ldr + i] + _cos * sdiag[i]
          r[(k-1) * ldr + i] = temp
        end
      end
    end

    sdiag[j] = r[(j -1) * ldr + j]
    r[(j-1) * ldr + j] = x[j]
  end

  # Solve the triangular system for z.  If the system is singular, then obtain a least-squares solution
  nsing = n
  for j=1:n
    if sdiag[j] == 0.0 && nsing ==n
      nsing = j - 1
    end

    if nsing < n
      W[j] = 0.0
    end
  end

  for j = nsing:-1:1
    _sum = 0.0
    for i = j+1:nsing # check this
      _sum += r[(j-1) * ldr + i] * W[i]
    end

    W[j] = (W[j] - _sum) / sdiag[j]
  end

  # permute the components of z back to the components of x
  for j=1:n
    x[pivot[j]] = W[j]
  end
end


function lm_lmpar!(n::Int, r::Matrix{Float64}, ldr::Int, pivot::Vector{Int}, diagonal::Vector{Float64}, qtb::Vector{Float64}, delta::Float64, par::Float64,
                  x::Vector{Float64}, sdiag::Vector{Float64}, aux::Vector{Float64}, xdi::Vector{Float64})

  # Compute and store in x the Gauss-Newton direction.  If the jacobian is rank-deficient, obtain a least-squares solution

  # println("here i am")

  nsing = n
  for j = 1:n
    aux[j] = qtb[j]
    if r[(j-1) * ldr + j] == 0.0 && nsing == n
      nsing = j - 1
    end

    if nsing < n
      aux[j] = 0.0
    end
  end

  for j = nsing:-1:1
    aux[j] = aux[j] / r[j + ldr * (j-1)]
    temp = aux[j]
    for i = 1:j - 1  ## CHECK THIS
      aux[i] -= r[(j-1) * ldr + i] * temp
    end
  end

  for j = 1:n
    x[pivot[j]] = aux[j]
  end

  # Initialize the iteration counter, evaluate the function at the origin,
  # and test for acceptance of the Gauss-Newton direction
  xdi = diagonal .* x

  dxnorm = vecnorm(xdi)
  fp = dxnorm - delta
  if fp <= 0.1 * delta
    par = 0.0
    return par
  end

  # if the jacobian is not rank deficient, the newton step provides a lower bound, parl, for the zero of the function.  Otherwise,
  # set this bound to zero

  parl = 0.0
  if nsing >= n
    for j =1:n
      aux[j] = diagonal[pivot[j]] * xdi[pivot[j]] / dxnorm
    end

    for j = 1:n
      _sum = 0.0
      for i = 1:j - 1 # CHECK THIS
        _sum += r[(j-1) * ldr + i] * aux[i]
      end

      aux[j] = (aux[j] - _sum) / r[j + ldr * (j-1)]
    end

    temp = vecnorm(aux)
    parl = fp / delta / temp / temp
  end

  # Calculate the upper bound, paru, for the zero of the function
  for j = 1:n
    _sum = 0.0
    for i = 1:j
      _sum += r[(j-1) * ldr + i] * qtb[i]
    end
    aux[j] = _sum / diagonal[pivot[j]]
  end

  # if the input par lies outside the interval (parl, paru), set par to the closer endpoint

  gnorm = vecnorm(aux)
  paru = gnorm / delta
  if paru == 0.0
    paru = DWARF / min(delta, 0.1)
  end

  # if the input par lies outside of the interval (parl, paru), set par to the closer endpoint

  par = max(par, parl)
  par = min(par, paru)
  if par == 0.0
    par = gnorm / dxnorm
  end

  # iterate
  for iter=1:10
    # evaluate the function at the current value of par
    if par == 0.0
      par = max(DWARF, 0.001 * paru)
    end

    temp = sqrt(par)
    aux = diagonal * temp

    lm_qrsolv!(n, r, ldr, pivot, aux, qtb, x, sdiag, xdi)

    xdi = diagonal .* x
    dxnorm = vecnorm(xdi)
    fp_old = fp
    fp = dxnorm - delta

    if abs(fp) <= 0.1 * delta || (parl == 0.0 && fp <= fp_old && fp_old < 0.0) || iter == 10
      break
    end

    # Compute the Newton correctly
    for j = 1:n
      aux[j] = diagonal[pivot[j]] * xdi[pivot[j]] / dxnorm
    end

    for j=1:n
      aux[j] = aux[j] / sdiag[j]
      for i=j+1:n
        aux[i] -= r[(j-1) * ldr + i] * aux[j]
      end
    end

    temp = vecnorm(aux)
    parc = fp / delta / temp / temp

    # depending on the sign of the function, update parl or paru
    if fp > 0
      parl = max(parl, par)
    else
      paru = min(paru, par)
    end

    par = max(parl, par * parc)
  end
  return par
end

function lmdif2!(n::Int, m::Int, x::Vector{Float64}, mode::Int, factor_::Float64, info_::Int, epsfcn::Float64, ftol::Float64, xtol::Float64, gtol::Float64,
                maxIter::Int, fcn!::Function, jacFcn!::Function = jacFcnDefault!)
  converged = false
  x_converged = false
  g_converged = false
  iterCount = 0
  f_calls = 0
  j_calls = 0
  eps_ = sqrt(max(epsfcn, MACHEP))

  # arrays to create
  fvec = zeros(m)
  diagonal = ones(n)
  qtf = zeros(n)
  fjac = zeros(m, n)
  wa1 = zeros(n)
  wa2 = zeros(n)
  wa3 = zeros(n)
  wf = zeros(m)
  pivot = ones(Int, n)

  # eval the function at the starting point and calculate its norm
  fcn!(m, n, x, fvec)
  f_calls += 1

  fnorm = vecnorm(fvec)

  # initialize the l-m parameter
  par = 0.0
  iterCount += 1

  while ~converged
    # Calculate the jacobian
    jacFcn!(m, n, x, fvec, fjac, eps_, wf, fcn!)
    j_calls += 1

    # println(fjac)

    qr = qrfact!(fjac, Val{true})
    # println(qr)
    # println("=======")
    # println(qr[:Q])
    # println("=======")
    # println(qr[:R])

    q = full(qr[:Q])
    r = qr[:R]
    pivot = qr[:p]
    wa1 = diag(r) # rdiag
    wa2 = r'[:, 1] # first of each column

    # println("fvec: ", fvec)
    #
    # println("next: ", q' * fvec)

    # Form q' * fvec and store in qtf
    #fjac = q
    wf = fjac * fvec
    qtf[1:n] = q' * fvec

    # println("wa1: ", wa1)
    # println("wa2: ", wa2)
    # println("qtf: ", qtf)
    # println("wf: ", wf)
    # println("fjac: ", fjac)
    #
    # error("break")

    # Compute norm of scaled gradient and detect degeneracy
    gnorm = 0.0
    for j = 1:n
      if wa2[pivot[j]] == 0.0
        continue
      end
      _sum = 0.0
      for i = 1:j
        _sum += fjac[i, j] * qtf[i]
      end

      gnorm = max(gnorm, abs(_sum / wa2[pivot[j]] / fnorm))
    end

    # terminate if gnorm <= gtol
    if gnorm <= gtol
      info_ = 4
      converged = true
      break
    end

    if iterCount == 1
      diagonal = copy(wa2)
      wa3 = diagonal .* x
      xnorm = vecnorm(wa3)

      delta = factor_ * xnorm
    else
      diagonal = max(diagonal, wa2)
    end

    inner_success = false

    while ~inner_success
      par = lm_lmpar!(n, fjac, m, pivot, diagonal, qtf, delta, par, wa1, wa2, wf, wa3)

      pnorm = vecnorm(wa3)

      temp2 = par * (pnorm / fnorm)^2

      for j = 1:n
        wa3[j] = 0.0
        for i = 1:j
          wa3[i] -= fjac[i, j] * wa1[pivot[j]]
        end
      end

      temp1 = (vecnorm(wa3) / fnorm)^2

      prered = temp1 + 2.0 * temp2
      dirder = -temp1 + temp2

      # at first call, adjust the initial step bound
      if iterCount == 1 && pnorm < delta
        delta = pnorm
      end

      # Evaluate the function at x + p
      wa2 = x - wa1
      fcn!(m, n, wa2, wf)
      f_calls += 1

      fnorm1 = vecnorm(wf)

      # Evaluate the scaled reduction
      actred = 1.0 - (fnorm1 / fnorm)^2

      # ratio of actual to predicted reduction
      ratio = prered > 0.0 ? actred / prered : 0.0

      # Update the step bound
      if ratio <= 0.25
        if actred >= 0.0
          temp = 0.5
        elseif actred > -99.0
          temp = max(dirder / (2.0 * dirder + actred), 0.1)
        else
          temp = 0.1
        end

        delta = temp * min(delta, pnorm / 0.1)
        par /= temp
      elseif ratio >= 0.75
        delta = 2.0 * pnorm
        par *= 0.5
      elseif par == 0.0
        delta = 2.0 * pnorm
      end

      # on success, update solution and test for convergence
      inner_success = ratio >= 1e-4

      if inner_success
        x = copy(wa2)
        wa2 = diagonal .* x
        fvec = copy(wf)

        xnorm = vecnorm(wa2)
        fnorm = fnorm1

        # Convergence tests
        info_ = 0

        if fnorm <= DWARF
          info_= 10 # success: sum of squares is almost 0
          converged = true
          break
        end

        if abs(actred) <= ftol && prered <= ftol && ratio <= 2.0
          info_ = 1 # success: x is almost stable
        end

        if delta <= xtol * xnorm
          info_ += 2 # success: sum of squares is almost stable
        end

        if info_ != 0
          converged = true
          break
        end

        # Tests for termination and stringent tolerances
        if abs(actred) <= MACHEP && prered <= MACHEP && ratio <= 2.0
          info_ = 6
          converged = true
          break
        end

        if delta <= MACHEP * xnorm
          info_ = 7
          converged = true
          break
        end

        if gnorm <= MACHEP
          info_ = 8
          converged = true
          break
        end
      end # convergence tests
    end # inner loop
    iterCount += 1
  end # outer loop
  return par, info_, f_calls, j_calls
end
