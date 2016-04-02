abstract EigenVectorCalculation
abstract ShiftStrategy

type WithEigenVector <: EigenVectorCalculation end
type WithoutEigenVector <: EigenVectorCalculation end
type OnlyFirstRowEigenVector <: EigenVectorCalculation end

type NoShift <: ShiftStrategy end
type Overrelaxation <: ShiftStrategy end
type CloseEigenValue <: ShiftStrategy end

immutable TqrEigenDecomposition
  iter::Int
  d::Vector{Float64}
  ev::Matrix{Float64}
end

function TqrEigenDecomposition(dg::Vector{Float64}, subVec::Vector{Float64}, calc::EigenVectorCalculation, strategy::ShiftStrategy)
  iter = 0
  d = dg
  ev_cols = get_ev_cols(calc)
  n = length(d)
  ev = zeros(n, ev_cols)
  e = copy(subVec)
  # insert!(e, 1, 0.0)
  unshift!(e, 0.0)

  for i = 1:ev_cols
    ev[i, i] = 1.0
  end

  @inbounds for k = n:-1:2
    while ~off_diag_is_zero(k, e, d)
      l = k
      l -= 1
      while l > 1 && ~off_diag_is_zero(l, e, d)
        l -= 1
      end

      iter += 1

      q = d[l]

      t1 = sqrt(0.25 * (d[k] * d[k] + d[k - 1] * d[k - 1]) - 0.5 * d[k - 1] * d[k] + e[k] * e[k])
      t2 = 0.5 * (d[k] + d[k - 1])

      lambda = abs(t2 + t1 - d[k]) < abs(t2 - t1 - d[k]) ? t2 + t1 : t2 - t1

      q -= (k == n ? 1.25 : 1.0) * lambda

      # QR transformation
      sine = 1.0
      cosine = 1.0
      u = 0.0

      recoverUnderflow = false

      @inbounds for i = l+1:k
        if recoverUnderflow
          break
        end

        h = cosine * e[i]
        p = sine * e[i]

        e[i - 1] = sqrt(p * p + q * q)
        if e[i - 1] != 0.0
          sine = p / e[i - 1]
          cosine = q / e[i - 1]

          g = d[i - 1] - u
          t = (d[i] - g) * sine + 2.0 * cosine * h

          u = sine * t
          d[i - 1] = g + u
          q = cosine * t - h

          @inbounds for j = 1:ev_cols
            tmp = ev[i - 1, j]
            ev[i - 1, j] = sine * ev[i, j] + cosine * tmp
            ev[i, j] = cosine * ev[i, j] - sine * tmp
          end
        else
          # recover from underflow
          d[i - 1] -= u
          e[l] = 0.0
          recoverUnderflow = true
        end
      end

      if ~recoverUnderflow
        d[k] -= u
        e[k] = q
        e[l] = 0.0
      end
    end
  end

  # sort (eigenvalues, eigenvectors)
  ev = ev[sortperm(d, rev=true), :]
  sort!(d, rev=true)
  @simd for i = 1:n
    sgn = 1.0
    @inbounds if ev_cols > 0 && ev[i, 1] < 0.0
      sgn = -1.0
    end

    @inbounds for j = 1:ev_cols
      ev[i, j] = sgn * ev[i, j]
    end
  end

  return TqrEigenDecomposition(iter, d, ev)
end

get_ev_cols(::OnlyFirstRowEigenVector) = 1

off_diag_is_zero(k::Int, e::Vector{Float64}, d::Vector{Float64}) = abs(d[k - 1]) + abs(d[k]) == abs(d[k - 1]) + abs(d[k]) + abs(e[k])
