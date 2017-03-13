immutable SVD
  U::Matrix{Float64}
  V::Matrix{Float64}
  s::Vector{Float64}
  transp::Bool
end

function svd_singular_val_creation!(::Type{Val{1}}, e::Vector{Float64}, s::Vector{Float64}, ::Matrix{Float64}, V::Matrix{Float64},
                                    p::Int, ::Int, iter::Int, k::Int, m::Int, n::Int)
  # Deflate negligible s(p)
  f = e[p-1]
  e[p-1] = 0.0
  @inbounds @simd for j = p-1:-1:k
    t = hypot(s[j],f)
    cs = s[j]/t
    sn = f/t
    s[j] = t
    if j != k
      f = -sn*e[j-1]
      e[j-1] = cs*e[j-1]
    end
    @inbounds for i = 1:n
      t = cs*V[i,j] + sn*V[i,p]
      V[i,p] = -sn*V[i,j] + cs*V[i,p]
      V[i,j] = t
    end
  end

  return iter, p, k
end

function svd_singular_val_creation!(::Type{Val{2}}, e::Vector{Float64}, s::Vector{Float64}, U::Matrix{Float64}, ::Matrix{Float64},
                                    p::Int, ::Int, iter::Int, k::Int, m::Int, n::Int)
  # Split at negligible s(k)
  f = e[k-1]
  e[k-1] = 0.0
  @inbounds @simd for j = k:p
    t = hypot(s[j],f)
    cs = s[j]/t
    sn = f/t
    s[j] = t
    f = -sn*e[j]
    e[j] = cs*e[j]
    @inbounds for i = 1:m
      t = cs*U[i,j] + sn*U[i,k-1]
      U[i,k-1] = -sn*U[i,j] + cs*U[i,k-1]
      U[i,j] = t
    end
  end
  return iter, p, k
end

function svd_singular_val_creation!(::Type{Val{3}}, e::Vector{Float64}, s::Vector{Float64}, U::Matrix{Float64}, V::Matrix{Float64},
                                    p::Int, ::Int, iter::Int, k::Int, m::Int, n::Int)
  # Perform one qr step

  # Calculate the shift
  scale_ = max(
            max(
              max(
                max(abs(s[p]), abs(s[p-1])),
                  abs(e[p-1])),
                abs(s[k])),
              abs(e[k]))

  sp = s[p]/scale_
  spm1 = s[p-1]/scale_
  epm1 = e[p-1]/scale_
  sk = s[k]/scale_
  ek = e[k]/scale_
  b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0
  c = (sp*epm1)*(sp*epm1)
  shift = 0.0
  if b != 0.0 || c != 0.0
    shift = sqrt(b*b + c)
    if b < 0.0
      shift = -shift
    end
    shift = c/(b + shift)
  end
  f = (sk + sp)*(sk - sp) + shift
  g = sk*ek

  # Chase zeros
  @inbounds @simd for j = k:p-1
    t = hypot(f,g)
    cs = f/t
    sn = g/t
    if j != k
      e[j-1] = t
    end
    f = cs*s[j] + sn*e[j]
    e[j] = cs*e[j] - sn*s[j]
    g = sn*s[j+1]
    s[j+1] = cs*s[j+1]
    @inbounds for i = 1:n
      t = cs*V[i,j] + sn*V[i,j+1]
      V[i,j+1] = -sn*V[i,j] + cs*V[i,j+1]
      V[i,j] = t
    end
    t = hypot(f,g)
    cs = f/t
    sn = g/t
    s[j] = t
    f = cs*e[j] + sn*s[j+1]
    s[j+1] = -sn*e[j] + cs*s[j+1]
    g = sn*e[j+1]
    e[j+1] = cs*e[j+1]
    if j < m
      @inbounds for i = 1:m
        t = cs*U[i,j] + sn*U[i,j+1];
        U[i,j+1] = -sn*U[i,j] + cs*U[i,j+1]
        U[i,j] = t
      end
    end
  end
  e[p-1] = f
  iter = iter + 1

  return iter, p, k
end

function svd_singular_val_creation!(::Type{Val{4}}, e::Vector{Float64}, s::Vector{Float64}, U::Matrix{Float64},
                                    V::Matrix{Float64}, p::Int, pp::Int, iter::Int, k::Int, m::Int, n::Int)
  # Make the singular values positive
  if s[k] <= 0.0
    s[k] = s[k] < 0.0 ? -s[k] : 0.0
    @simd for i = 1:pp+1
      @inbounds V[i,k] = -V[i,k]
    end
  end

  # Order the singular values
  @inbounds while k <= pp
    if s[k] >= s[k+1]
      break
    end
    s[k], s[k+1] = s[k+1], s[k]
    if k < n
      @inbounds for i = 1:n
        V[i, k], V[i, k+1] = V[i, k+1], V[i, k]
      end
    end
    if k < m
      @inbounds for i = 1:m
        U[i, k], U[i, k+1] = U[i, k+1], U[i, k]
      end
    end
    k += 1
  end

  iter = 0
  p -= 1

  return iter, p, k
end

function SVD(M::Matrix{Float64})
  m, n = size(M)

  if m >= n
    A = copy(M)
    transp = false
  else
    A = M'
    transp = true
  end

  s = Vector{Float64}(n)
  U = zeros(m, n)
  V = Matrix{Float64}(n, n)
  e = Vector{Float64}(n)
  work = Vector{Float64}(m)
  i = j = k = 0

  # reduce A to its bidiagonal form, storing the diagonal elements
  # in s, and super-diagonal elements in e
  nct = min(m-1, n)
  nrt = max(0, n-2)

  @inbounds @simd for k = 1:max(nct, nrt)
    if k <= nct
      # Compute the transformation for the k-th column and
      # place the k-th diagonal in s[k].
      # Compute 2-norm of k-th column without under/overflow.
      s[k] = 0.0
      @inbounds for i = k:m
        s[k] = hypot(s[k],A[i,k]);
      end
      if s[k] != 0.0
        if A[k,k] < 0.0
          s[k] = -s[k];
        end
        for i = k:m
          A[i,k] /= s[k]
        end
        A[k,k] += 1.0;
      end
      s[k] = -s[k]
    end
    @inbounds for j = k+1:n
      if k <= nct && s[k] != 0.0
        # apply transformation
        t = 0.0
        @inbounds for i = k:m
          t += A[i,k]*A[i,j]
        end
        t = -t/A[k,k]
        @inbounds for i = k:m
          A[i,j] += t*A[i,k]
        end
      end
      # Place the k-th row of A into e for the
      # subsequent calculation of the row transformation.
      e[j] = A[k, j]
    end
    if k < nct
      # Place the transformation in U for subsequent back
      # multiplication.
      @inbounds for i = k:m
        U[i,k] = A[i,k]
      end
    end
    if k <= nrt
      # Compute the k-th row transformation and place the
      # k-th super-diagonal in e[k].
      # Compute 2-norm without under/overflow.
      e[k] = 0.0
      @inbounds for i = k+1:n
        e[k] = hypot(e[k],e[i])
      end
      if e[k] != 0.0
        if e[k+1] < 0.0
          e[k] = -e[k]
        end
        @inbounds for i = k+1:n
          e[i] /= e[k]
        end
        e[k+1] += 1.0
      end
      e[k] = -e[k]
      if k+1 <= m && e[k] != 0.0

        # Apply the transformation.
        for i = k+1:m
          work[i] = 0.0
        end
        @inbounds for j = k+1:n
          @inbounds for i = k+1:m
            work[i] += e[j]*A[i,j]
          end
        end
        @inbounds for j = k+1:n
          t = -e[j]/e[k+1]
          @inbounds for i = k+1:m
            A[i,j] += t*work[i]
          end
        end
      end

      # Place the transformation in V for subsequent
      # back multiplication.
      @inbounds for i = k+1:n
        V[i,k] = e[i]
      end
    end
  end

  # Set up the final bidiagonal matrix or order n.

  if nct < n
    s[nct+1] = A[nct+1, nct+1]
  end
  if nrt+1 < n
    e[nrt+1] = A[nrt+1, n]
  end
  e[n] = 0.0

  # Generate U
  @inbounds @simd for j = nct+1:n
    for i = 1:m
      U[i,j] = 0.0
    end
    U[j,j] = 1.0
  end

  @inbounds @simd for k = nct:-1:1
    if s[k] != 0.0
      @inbounds for j = k+1:n
        t = 0
        @inbounds for i = k:m
          t += U[i,k]*U[i,j]
        end
        t = -t/U[k,k]
        @inbounds for i = k:m
          U[i,j] += t*U[i,k]
        end
      end
      @inbounds for i = k:m
        U[i,k] = -U[i,k]
      end
      U[k,k] = 1.0 + U[k,k]
      @inbounds for i = 1:k-1
        U[i,k] = 0.0
      end
    else
      @inbounds for i = 1:m
        U[i,k] = 0.0
      end
      U[k,k] = 1.0
    end
  end

  # Generate V
  @inbounds @simd for k = n:-1:1
    if k <= nrt && e[k] != 0.0
      @inbounds for j = k+1:n
        t = 0.0
        @inbounds for i = k+1:n
          t += V[i,k]*V[i,j]
        end
        t = -t/V[k+1,k]
        @inbounds for i = k+1:n
          V[i,j] += t*V[i,k]
        end
      end
    end
    @inbounds for i = 1:n
      V[i,k] = 0.0
    end
    V[k,k] = 1.0
  end

  # Main iteration loop for singular values
  p = n
  pp = p-1
  iter = 0

  @inbounds while p > 0
    # Here is where a test for too many iterations would go.

    # This section of the program inspects for
    # negligible elements in the s and e arrays.  On
    # completion the variables kase and k are set as follows.

    # kase = 1     if s(p) and e[k-1] are negligible and k<p
    # kase = 2     if s(k) is negligible and k<p
    # kase = 3     if e[k-1] is negligible, k<p, and
    #              s(k), ..., s(p) are not negligible (qr step).
    # kase = 4     if e(p-1) is negligible (convergence).

    for k = p -1:-1:0
      if k == 0
        break
      end

      if abs(e[k]) <= EPS_VAL*(abs(s[k]) + abs(s[k+1]))
        e[k] = 0.0
        break
      end
    end

    if k == p-1
        kase = 4
    else
      ks = 0
      @inbounds for ks = p:-1:k
        if ks == k
          break
        end

        t = (ks != p+1 ? abs(e[ks]) : 0.0) + (ks != k+1 ? abs(e[ks-1]) : 0.0)
        if abs(s[ks]) <= EPS_VAL * t
          s[ks] = 0.0;
          break
        end
      end
      if ks == k
        kase = 3
      elseif ks == p
        kase = 1
      else
        kase = 2
        k = ks
      end
    end
    k += 1

    # perform the task indicated by kase
    iter, p, k = svd_singular_val_creation!(Val{kase}, e, s, U, V, p, pp, iter, k, m, n)

  end
  return SVD(U, V, s, transp)
end
