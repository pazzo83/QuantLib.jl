using LinearAlgebra

mutable struct TridiagonalOperator
  diagonal::Vector{Float64}
  lowerDiagonal::Vector{Float64}
  upperDiagonal::Vector{Float64}
  temp::Vector{Float64}
  n::Int
end

TridiagonalOperator(n::Int) = TridiagonalOperator(zeros(n), zeros(n - 1), zeros(n - 1), zeros(n), n)
TridiagonalOperator() = TridiagonalOperator(Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), Vector{Float64}(), 0)
TridiagIdentity(n::Int) = TridiagonalOperator(ones(n), zeros(n - 1), zeros(n - 1), zeros(n), n)

function set_first_row!(L::TridiagonalOperator, valB::Float64, valC::Float64)
  L.diagonal[1] = valB
  L.upperDiagonal[1] = valC

  return L
end

function set_mid_row!(L::TridiagonalOperator, i::Int, valA::Float64, valB::Float64, valC::Float64)
  i >= 2 && i <= L.n - 1 || error("out of range for Tridiagonal System: set mid row")

  L.lowerDiagonal[i - 1] = valA
  L.diagonal[i] = valB
  L.upperDiagonal[i] = valC

  return L
end

function set_last_row!(L::TridiagonalOperator, valA::Float64, valB::Float64)
  L.lowerDiagonal[L.n - 1] = valA
  L.diagonal[L.n] = valB

  return L
end

function solve_for(L::TridiagonalOperator, rhs::Vector{Float64})
  # create tridiagonal
  tri = Tridiagonal(L.lowerDiagonal, L.diagonal, L.upperDiagonal)

  # lu factorization
  Alu = lu(tri)

  return Alu \ rhs
end

function solve_for!(L::TridiagonalOperator, rhs::Vector{Float64}, result::Vector{Float64})
  bet = L.diagonal[1]
  result[1] = rhs[1] / bet

  @inbounds @simd for j = 2:L.n
    L.temp[j] = L.upperDiagonal[j - 1] / bet
    bet = L.diagonal[j] - L.lowerDiagonal[j - 1] * L.temp[j]
    result[j] = (rhs[j] - L.lowerDiagonal[j - 1] * result[j - 1]) / bet
  end

  @inbounds @simd for j = L.n-1:-1:2
    result[j] -= (L.temp[j + 1] * result[j + 1])
    # println(result[j])
  end
  # println(result)
  # error("DIE")
  result[1] -= L.temp[2] * result[2]
  return L, result
end

function apply_to(L::TridiagonalOperator, v::Vector{Float64})
  L.n == length(v) || error("v of incorrect length")

  res = L.diagonal .* v

  res[1] = L.upperDiagonal[1] * v[2]
  @simd for j = 2:L.n - 1
    @inbounds res[j] += L.lowerDiagonal[j-1] * v[j-1] + L.upperDiagonal[j] * v[j+1]
  end

  res[end] += L.lowerDiagonal[end-1] * v[end-1]

  return res
end

# basic operators
import Base.+, Base.-, Base.*, Base./, Base.copy

function -(D::TridiagonalOperator)
  low = -D.lowerDiagonal
  mid = -D.diagonal
  high = -D.upperDiagonal

  return TridiagonalOperator(mid, low, high, D.temp, D.n)
end

function +(D1::TridiagonalOperator, D2::TridiagonalOperator)
  D1.n == D2.n || error("Dimension mismatch")
  low = D1.lowerDiagonal + D2.lowerDiagonal
  mid = D1.diagonal + D2.diagonal
  high = D1.upperDiagonal + D2.upperDiagonal

  return TridiagonalOperator(mid, low, high, D1.temp, D1.n)
end

function -(D1::TridiagonalOperator, D2::TridiagonalOperator)
  D1.n == D2.n || error("Dimension mismatch")
  low = D1.lowerDiagonal - D2.lowerDiagonal
  mid = D1.diagonal - D2.diagonal
  high = D1.upperDiagonal - D2.upperDiagonal

  return TridiagonalOperator(mid, low, high, D1.temp, D1.n)
end

function *(a::Number, D::TridiagonalOperator)
  low = D.lowerDiagonal * a
  mid = D.diagonal * a
  high = D.upperDiagonal * a

  return TridiagonalOperator(mid, low, high, D.temp, D.n)
end

function *(D::TridiagonalOperator, a::Number)
  low = D.lowerDiagonal * a
  mid = D.diagonal * a
  high = D.upperDiagonal * a

  return TridiagonalOperator(mid, low, high, D.temp, D.n)
end

function /(D::TridiagonalOperator, a::Number)
  low = D.lowerDiagonal / a
  mid = D.diagonal / a
  high = D.upperDiagonal / a

  return TridiagonalOperator(mid, low, high, D.temp, D.n)
end

copy(L::TridiagonalOperator) = TridiagonalOperator(L.diagonal, L.lowerDiagonal, L.upperDiagonal, L.temp, L.n)
