const MACHEP = 1.2e-16
const DWARF = 1.0e-38

function enorm(n::Int, x::Vector{Float64})
  rdwarf = 3.834e-20
  rgiant = 1.304e19

  s1 = s2 = s3 = x1max = x3max = 0.0

  floatn = float(n)
  agiant = rgiant / floatn

  for i = 1:n
    @inbounds xabs = abs(x[i])
    if xabs > rdwarf && xabs < agiant
      # sum for intermediate components
      s2 += xabs * xabs
      continue
    end

    if xabs > rdwarf
      # sum for larger components
      if xabs > x1max
        temp = x1max / xabs
        s1 = 1.0 + s1 * temp * temp
        x1max = xabs
      else
        temp = xabs / x1max
        s1 += temp * temp
      end

      continue
    end

    # sum for smaller components
    if xabs > x3max
      temp = x3max / xabs
      s3 = 1.0 + s3 * temp * temp
      x3max = xabs
    else
      if xabs != 0.0
        temp = xabs / x3max
        s3 += temp * temp
      end
    end
  end

  # Calculation of norm
  if s1 != 0.0
    temp = s1 + (s2 / x1max) / x1max
    return x1max * sqrt(temp)
  end

  if s2 != 0.0
    if s2 >= x3max
      temp = s2 * (1.0 + (x3max / s2) * (x3max * s3))
    else
      temp = x3max * ((s2 / x3max) + (x3max * s3))
    end
    ans = sqrt(temp)
  else
    ans = x3max * sqrt(s3)
  end

  return ans
end

# Forward Difference Approximation #
function fdjac2!(m::Int, n::Int, x::Vector{Float64}, fvec::Vector{Float64}, fjac::Matrix{Float64}, ::Int, iflag::Int, epsfcn::Float64, wa::Vector{Float64}, fcn!::Function)
  # returns:
  # fjac :  is an output m by n array which contains the
  #         approximation to the jacobian matrix evaluated at x.
  temp  = max(epsfcn, MACHEP)
  eps_ = sqrt(temp)
  ij = 1
  @simd for j = 1:n
    @inbounds temp = x[j]
    h = eps_ * abs(temp)
    if h == 0.0
      h = eps_
    end

    @inbounds x[j] = temp + h
    fcn!(m, n, x, wa)
    if iflag < 0
      return fjac
    end
    @inbounds x[j] = temp
    for i = 1:m
      @inbounds fjac[ij] = (wa[i] - fvec[i]) / h
      ij += 1
    end
  end

  return fjac
end

function qrfac!(m::Int, n::Int, a::Matrix{Float64}, ::Int, pivot::Int, ipvt::Vector{Int}, ::Int, rdiag::Vector{Float64}, acnorm::Vector{Float64}, wa::Vector{Float64})
  # returns:
  # a :     is an m by n array. on input a contains the matrix for
  #         which the qr factorization is to be computed. on output
  #         the strict upper trapezoidal part of a contains the strict
  #         upper trapezoidal part of r, and the lower trapezoidal
  #         part of a contains a factored form of q (the non-trivial
  #         elements of the u vectors described above).
  # ipvt:   is an integer output array of length lipvt. ipvt
  #         defines the permutation matrix p such that a*p = q*r.
  #         column j of p is column ipvt(j) of the identity matrix.
  #         if pivot is false, ipvt is not referenced.
  # rdiag:  is an output array of length n which contains the
  #         diagonal elements of r.
  # acnorm: is an output array of length n which contains the
  #         norms of the corresponding columns of the input matrix a.
  #         if this information is not needed, then acnorm can coincide
  #         with rdiag.
  # Compute the initial column norms and initialize several arrays
  ij = 1
  @simd for j = 1:n
    @inbounds acnorm[j] = enorm(m, a[:, ij])
    @inbounds rdiag[j] = acnorm[j]
    @inbounds wa[j] = rdiag[j]
    if pivot != 0
      @inbounds ipvt[j] = j
    end
    ij += 1
  end

  # Reduce a to r with householder transformations
  minmn = min(m, n)
  @simd for j = 1:minmn
    if pivot != 0
      # bring the column of the largest norm into the pivot position
      kmax = j
      for k = j:n
        @inbounds if rdiag[k] > rdiag[kmax]
          kmax = k
        end
      end

      if kmax != j
        ij = (m * (j - 1)) + 1
        jj = (m * (kmax - 1)) + 1

        for i = 1:m
          @inbounds temp = a[ij]
          @inbounds a[ij] = a[jj]
          @inbounds a[jj] = temp
          ij += 1
          jj += 1
        end

        @inbounds rdiag[kmax] = rdiag[j]
        @inbounds wa[kmax] = wa[j]
        @inbounds k = ipvt[j]
        @inbounds ipvt[j] = ipvt[kmax]
        @inbounds ipvt[kmax] = k
      end
    end

    # Compute the householder transformation to reduce the j-th column of a to
    # a multiple of the j-th unit vector
    jj = j + m * (j - 1)
    # jj_end = jj + (m - j)
    @inbounds ajnorm = enorm(m - (j - 1), a[jj:end])
    if ajnorm != 0.0
      @inbounds if a[jj] < 0.0
        ajnorm = -ajnorm
      end
      ij = jj

      for i = j:m
        @inbounds a[ij] /= ajnorm
        ij += 1

      end

      @inbounds a[jj] += 1.0

      # Apply the transformation to the remaining columns
      # and update the norms
      jp1 = j + 1
      if jp1 <= n
        for k = jp1:n
          sum_ = 0.0
          ij = j + m * (k - 1)
          jj = j + m * (j - 1)

          for i = j:m
            @inbounds sum_ += a[jj] * a[ij]
            ij += 1
            jj += 1
          end

          @inbounds temp = sum_ / a[j + m * (j - 1)]
          ij = j + m * (k - 1)
          jj = j + m * (j - 1)
          for i = j:m
            @inbounds a[ij] -= temp * a[jj]
            ij += 1
            jj += 1
          end

          if pivot != 0 && rdiag[k] != 0.0
            @inbounds temp = a[j + m * (k - 1)]  / rdiag[k]
            temp = max(0.0, 1.0 - temp * temp)
            @inbounds rdiag[k] *= sqrt(temp)
            @inbounds temp = rdiag[k]/wa[k]

            if 0.05 * temp * temp <= MACHEP
              a_start = jp1 + m * (k - 1)
              @inbounds rdiag[k] = enorm(m - j, a[a_start:end])
              @inbounds wa[k] = rdiag[k]
            end
          end
        end
      end
    end

    @inbounds rdiag[j] = -ajnorm
  end

  return a, ipvt, rdiag, acnorm
end

function qrsolv!(n::Int, r::Matrix{Float64}, ldr::Int, ipvt::Vector{Int}, diag_::Vector{Float64}, qtb::Vector{Float64}, x::Vector{Float64},
                sdiag::Vector{Float64}, wa::Vector{Float64})
  # copy r and (q transpose) * b to preserve input and initialize s, in particular
  # save the diagonal elements of r in x
  kk = 1
  @simd for j = 1:n
    ij = kk
    ik = kk
    for i = j:n
      @inbounds r[ij] = r[ik]
      ij += 1
      ik += ldr
    end
    @inbounds x[j] = r[kk]
    @inbounds wa[j] = qtb[j]
    kk += ldr + 1
  end

  # eliminate the diagonal matrix d using a given rotation
  for j = 1:n
    @inbounds l = ipvt[j]
    @inbounds if diag_[l] != 0.0
      for k = j:n
        @inbounds sdiag[k] = 0.0
      end
      @inbounds sdiag[j] = diag_[l]

      # transformations to eliminate the row of d, modify only a single element
      # of (q transpose) * b beyond the first n, which is initially zero
      qtbpj = 0.0
      for k = j:n
        # determine a given rotation which elminates the appropriate element
        # in the current row of d
        @inbounds if sdiag[k] == 0.0
          continue
        end
        kk = k + ldr * (k - 1)
        @inbounds if abs(r[kk]) < abs(sdiag[k])
          @inbounds cotan_ = r[kk] / sdiag[k]
          sin_ = 0.5 / sqrt(0.25 + 0.25 * cotan_ * cotan_)
          cos_ = sin_ * cotan_
        else
          @inbounds tan_ = sdiag[k] / r[kk]
          cos_ = 0.5 / sqrt(0.25 + 0.25 * tan_ * tan_)
          sin_ = cos_ * tan_
        end

        # compute the modified diagonal element of r and the modified element of
        # (q transpose) * b, 0
        @inbounds r[kk] = cos_ * r[kk] + sin_* sdiag[k]
        @inbounds temp = cos_ * wa[k] + sin_ * qtbpj
        @inbounds qtbpj = -sin_ * wa[k] + cos_ * qtbpj
        @inbounds wa[k] = temp

        # accumulate the transformation of the row s
        kp1 = k + 1
        if n >= kp1
          ik = kk + 1
          for i = kp1:n
            @inbounds temp = cos_ * r[ik] + sin_ * sdiag[i]
            @inbounds sdiag[i] = -sin_ * r[ik] + cos_ * sdiag[i]
            @inbounds r[ik] = temp
            ik += 1
          end
        end
      end
    end

    # store the diagonal element of s and restore the corresponding diagonal
    # element of r
    kk = j + ldr * (j - 1)
    @inbounds sdiag[j] = r[kk]
    @inbounds r[kk] = x[j]
  end

  # Solve the triangular system for z, if the system is singular, then obtain
  # a least squares solution
  nsing = n
  @simd for j = 1:n
    @inbounds if sdiag[j] == 0.0 && nsing == n
      nsing = j - 1
    end

    if nsing < n
      @inbounds wa[j] = 0.0
    end
  end

  if nsing >= 1
    @simd for k = 1:nsing
      j = nsing - k + 1
      sum_ = 0.0
      jp1 = j + 1
      if nsing >= jp1
        ij = jp1 + ldr * (j - 1)
        for i = jp1:nsing
          @inbounds sum_ += r[ij] * wa[i]
          ij += 1
        end
      end
      @inbounds wa[j] = (wa[j] - sum_) / sdiag[j]
    end
  end
  # permute the components of z back to components of x
  @simd for j = 1:n
    # l = ipvt[j]
    @inbounds x[ipvt[j]] = wa[j]
  end

  return r, x, sdiag
end

function lmpar!(n::Int, r::Matrix{Float64}, ldr::Int, ipvt::Vector{Int}, diag_::Vector{Float64}, qtb::Vector{Float64}, delta::Float64, par::Float64,
              x::Vector{Float64}, sdiag::Vector{Float64}, wa1::Vector{Float64}, wa2::Vector{Float64})

  # Compute and store in x the gauss-newton direction.  If the jacobin is rank-deficient
  # obtain a least-squares solution

  nsing = n
  jj = 1
  @simd for j=1:n
    @inbounds wa1[j] = qtb[j]
    @inbounds if r[jj] == 0.0 && nsing == n
      nsing = j - 1
    end
    if nsing < n
      @inbounds wa1[j] = 0.0
    end

    jj += ldr + 1
  end

  if nsing >= 1
    @simd for k = 1:nsing
      j = nsing - k + 1
      @inbounds wa1[j] = wa1[j] / r[j + ldr * (j - 1)]
      @inbounds temp = wa1[j]
      jm1 = j - 1
      if jm1 > 0
        ij = ldr * (j - 1) + 1
        for i = 1:jm1 + 1
          @inbounds wa1[i] -= r[ij] * temp
          ij += 1
        end
      end
    end
  end

  @simd for j = 1:n
    @inbounds l = ipvt[j]
    @inbounds x[l] = wa1[j]
  end

  # initialize the iteration counter, evaluate the function at the origin
  # and test for acceptance of the gauss-newton direction
  iter = 0
  @simd for j = 1:n
    @inbounds wa2[j] = diag_[j] * x[j]
  end

  dxnorm = enorm(n, wa2)
  fp = dxnorm - delta

  # fp <= 0.1 * delta && @goto L220
  if fp > 0.1 * delta
    # if the jacobin is not rank deficient, the newton step provides a lower bound
    # parl, for the zero of the function.  Otherwise set this bound to zero
    parl = 0.0
    if nsing >= n
      @simd for j = 1:n
        @inbounds l = ipvt[j]
        @inbounds wa1[l] = diag_[l] * (wa2[l] / dxnorm)
      end
      jj = 1
      @simd for j = 1:n
        sum_ = 0.0
        jm1 = j - 1
        if jm1 >= 1
          ij = jj
          for i = 1:jm1 + 1
            @inbounds sum_ += r[ij] * wa1[i]
            ij += 1
          end
        end
        @inbounds wa1[j] = (wa1[j] - sum_) / r[j + ldr * (j - 1)]
        jj += ldr
      end
      temp = enorm(n, wa1)
      parl = ((fp / delta) / temp) / temp
    end

    jj = 1
    @simd for j = 1:n
      sum_ = 0.0
      ij = jj
      for i = 1:j
        @inbounds sum_ += r[ij] * qtb[i]
        ij += 1
      end
      @inbounds l = ipvt[j]
      @inbounds wa1[j] = sum_ / diag_[l]
      jj += ldr
    end

    gnorm = enorm(n, wa1)
    paru = gnorm / delta
    if paru == 0.0
      paru = DWARF / min(delta, 0.1)
    end

    # if the input par lies outside of the interval (parl, paru), set the par to
    # the closer endpoint
    par = max(par, parl)
    par = min(par, paru)

    if par == 0.0
      par = gnorm / dxnorm
    end
    continueiter = true
    while continueiter
      iter += 1

      # eval the function at the current value of par
      if par == 0.0
        par = max(DWARF, 0.001 * paru)
      end

      temp = sqrt(par)
      # for j = 1:n
      #   wa1[j] = temp * diag_[j]
      # end
      wa1 = temp * diag_

      qrsolv!(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2)

      @simd for j = 1:n
        @inbounds wa2[j] = diag_[j] * x[j]
      end

      dxnorm = enorm(n, wa2)
      temp = fp
      fp = dxnorm - delta

      # if the function is small enough, accept the current value of par.  also test
      # for the exceptional cases where parl is zero or the number of iterations
      # has reached 10
      # (abs(fp) <= 0.1 * delta || (parl == 0.0 && fp <= temp && temp < 0.0) || iter == 10) && @goto L220
      if (abs(fp) <= 0.1 * delta || (parl == 0.0 && fp <= temp && temp < 0.0) || iter == 10)
        continueiter = false
      else
        # compute the newton correction
        @simd for j = 1:n
          @inbounds l = ipvt[j]
          @inbounds wa1[j] = diag_[l] * (wa2[l] / dxnorm)
        end

        jj = 1
        @simd for j = 1:n
          @inbounds wa1[j] = wa1[j] / sdiag[j]
          @inbounds temp = wa1[j]
          jp1 = j + 1
          if jp1 <= n
            ij = jp1 + jj
            for i = jp1:n
              @inbounds wa1[i] -= r[ij] * temp
              ij += 1
            end
          end
          jj += ldr
        end
        temp = enorm(n, wa1)
        parc = ((fp / delta) / temp) / temp

        # depending on the sign of the function, update parl or paru
        if fp > 0.0
          parl = max(parl, par)
        end

        if fp < 0.0
          paru = min(paru, par)
        end

        # compute an improved estimate for par
        par = max(parl, par + parc)
      end
      # end of iteration
    end
  end
  # termination
  if iter == 0
    par = 0.0
  end

  return par
end

function lmdif!(m::Int, n::Int, x::Vector{Float64}, fvec::Vector{Float64}, ftol::Float64, xtol::Float64, gtol::Float64, maxFev::Int, epsfcn::Float64,
                diag_::Vector{Float64}, mode::Int, factor_::Float64, nprint::Int, info_::Int, nfev::Int, fjac::Matrix{Float64}, ldfjac::Int, ipvt::Vector{Int}, qtf::Vector{Float64},
                wa1::Vector{Float64}, wa2::Vector{Float64}, wa3::Vector{Float64}, wa4::Vector{Float64}, fcn!::Function)

  delta = 0.0
  xnorm = 0.0
  info_ = 0
  iflag = 0
  nfev = 0

  # checking for errors
  mainloop = true
  while mainloop
    mainloop = false # only iterate once
    (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0 || gtol < 0.0 || maxFev <= 0 || factor_ <= 0.0) && break

    if mode == 2
      for j = 1:n
        @inbounds if diag_[j] <= 0.0
          break
        end
      end
    end

    # Evaluate the function at its starting point and calculate its norm
    iflag = 1
    fcn!(m, n, x, fvec)
    nfev = 1
    if iflag < 0
      break
    end
    fnorm = enorm(m, fvec)

    # Initialize the levenberg-marquardt param and iteration counter
    par = 0.0
    iter = 1

    outer = true
    while outer
      # Calcualte the Jacobian matrix
      iflag = 2
      # if useJac
      #   jacFcn!(m, n, x, fjac, iflag)
      # else
      fdjac2!(m, n, x, fvec, fjac, ldfjac, iflag, epsfcn, wa4, fcn!)
      # end
      # fjac = jacFcn!(x)

      nfev += n

      if iflag < 0
        break
      end

      # if requested, call fcn to enable printing of iterates
      if nprint > 0
        iflag = 0
        if mod(iter - 1, nprint) == 0
          fcn!(m, n, x, fvec)
          if iflag < 0
            break
          end
        end
      end

      # compute the QR factorization of the Jacobian
      qrfac!(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3)

      # on the first iteration and if mode is 1, scale according to the norms
      # of the columns in the initial jacobin
      if iter == 1
        if mode != 2
          @simd for j = 1:n
            @inbounds diag_[j] = wa2[j]
            if wa2[j] == 0.0
              @inbounds diag_[j] = 1.0
            end
          end
        end

        # on the first iteration, calculate the norm of the scaled x and initialize
        # the step bound delta
        @simd for j = 1:n
          @inbounds wa3[j] = diag_[j] * x[j]
        end

        xnorm = enorm(n, wa3)
        delta = factor_ * xnorm
        if delta == 0.0
          delta = factor_
        end

      end

      # form (q transpose) * fvec and store the first n components in qtf
      @simd for i = 1:m
        @inbounds wa4[i] = fvec[i]
      end

      jj = 1
      @simd for j = 1:n
        @inbounds temp3 = fjac[jj]
        if temp3 != 0.0
          sum_ = 0.0
          ij = jj
          for i = j:m
            @inbounds sum_ += fjac[ij] * wa4[i]
            ij += 1
          end
          temp = -sum_ / temp3
          ij = jj
          for i = j:m
            @inbounds wa4[i] += fjac[ij] * temp
            ij += 1
          end
        end
        @inbounds fjac[jj] = wa1[j]
        jj += m + 1
        @inbounds qtf[j] = wa4[j]
      end


      # compute the norm of the scaled gradient
      gnorm = 0.0
      if fnorm != 0.0
        jj = 1
        @simd for j = 1:n
          @inbounds l = ipvt[j]
          @inbounds if wa2[l] != 0.0
            sum_ = 0.0
            ij = jj
            for i = 1:j
              @inbounds sum_ += fjac[ij] * (qtf[i] / fnorm)
              ij += 1
            end
            @inbounds gnorm = max(gnorm, abs(sum_ / wa2[l]))
          end
          jj += m
        end
      end

      # test for convergence of the gradient norm
      if gnorm <= gtol
        info_ = 4
      end

      if info_ != 0
        break
      end

      # rescale if necessary
      if mode != 2
        @simd for j = 1:n
          @inbounds diag_[j] = max(diag_[j], wa2[j])
        end
      end

      # Beginning of inner loop
      inner = true
      while inner

        # determine the levenberg-marquardt param
        par = lmpar!(n, fjac, ldfjac, ipvt, diag_, qtf, delta, par, wa1, wa2, wa3, wa4)

        # store the direction of p and x + p, calculate the norm of p
        @simd for j = 1:n
          @inbounds wa1[j] = -wa1[j]
          @inbounds wa2[j] = x[j] + wa1[j]
          @inbounds wa3[j] = diag_[j] * wa1[j]
        end

        pnorm = enorm(n, wa3)

        # on the initial iteration, adjust the initial step bound
        if iter == 1
          delta = min(delta, pnorm)
        end

        iflag = 1
        fcn!(m, n, wa2, wa4)
        nfev += 1

        if iflag < 0
          outer = false
          break
        end

        fnorm1 = enorm(m, wa4)
        # println(wa4)

        # Compute the scaled actual reduction
        actred = -1.0

        # println("fnorm: ", fnorm)
        # println("fnorm1: ", fnorm1)

        if (0.1 * fnorm1) < fnorm
          temp = fnorm1 / fnorm
          actred = 1.0 - temp * temp
        end

        # Compute the scaled predicted reduction and the scaled directional derivative
        jj = 1
        @simd for j = 1:n
          @inbounds wa3[j] = 0.0
          @inbounds l = ipvt[j]
          @inbounds temp = wa1[l]
          ij = jj
          for i = 1:j
            @inbounds wa3[i] += fjac[ij] * temp
            ij += 1
          end

          jj += m
        end

        temp1 = enorm(n, wa3) / fnorm
        temp2 = (sqrt(par) * pnorm) / fnorm
        prered = temp1 * temp1 + (temp2 * temp2) / 0.5
        dirder = -(temp1 * temp1 + temp2 * temp2)

        # Compute the ratio of the actual to the predicted reduction
        ratio = 0.0
        if prered != 0.0
          ratio = actred / prered
        end

        if ratio <= 0.25
          if actred >= 0.0
            temp = 0.5
          else
            temp = 0.5 * dirder / (dirder + 0.5 * actred)
          end

          if ((0.1 * fnorm1) >= fnorm) || temp < 0.1
            temp = 0.1
          end

          delta = temp * min(delta, pnorm/0.1)
          par = par / temp
        else
          if par == 0.0  || ratio >= 0.75
            delta = pnorm / 0.5
            par = 0.5 * par
          end
        end
        # test for successful iteration
        if ratio >= 0.0001
          @simd for j = 1:n
            @inbounds x[j] = wa2[j]
            @inbounds wa2[j] = diag_[j] * x[j]
          end
          @simd for i = 1:m
            @inbounds fvec[i] = wa4[i]
          end

          xnorm = enorm(n, wa2)
          fnorm = fnorm1
          iter += 1
        end

        # Tests for convergence
        if abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1.0
          info_ = 1
        end

        if delta <= xtol * xnorm
          info_ = 2
        end

        if abs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1.0 && info_ == 2
          info_ = 3
        end

        if info_ != 0
          outer = false
          break
        end

        # Tests for termination and stringent tolerances
        if nfev >= maxFev
          info_ = 5
        end

        if abs(actred) <= MACHEP && prered <= MACHEP && 0.5 * ratio <= 1.0
          info_ = 6
        end

        if delta <= MACHEP * xnorm
          info_ = 7
        end

        if gnorm <= MACHEP
          info_ = 8
        end

        if info_ != 0
          outer = false
          break
        end

        # end of inner loop, repeat if iteration unsuccessful
        if ratio >= 0.0001
          inner = false
        end
      end
    end
    # end of outer loop
    # @goto L30
  end

  # @label L300
  # Termination either normal or user imposed
  if iflag < 0
    info_ = iflag
  end
  iflag = 0
  if nprint > 0
    fcn!(m, n, x, fvec)
  end

  return x, fvec, info_, nfev, fjac, ipvt, qtf
end
