type TridiagonalOperator{I <: Integer}
  diagonal::Vector{Float64}
  lowerDiagonal::Vector{Float64}
  upperDiagonal::Vector{Float64}
  temp::Vector{Float64}
  n::I
end

TridiagonalOperator(n::Int) = TridiagonalOperator(zeros(n), zeros(n - 1), zeros(n - 1), zeros(n), n)

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
  Alu = lufact(tri)

  return Alu \ rhs
end

function solve_for!(L::TridiagonalOperator, rhs::Vector{Float64}, result::Vector{Float64})
  bet = L.diagonal[1]
  result[1] = rhs[1] / bet

  for j = 2:L.n
    L.temp[j] = L.upperDiagonal[j - 1] / bet
    bet = L.diagonal[j] - L.lowerDiagonal[j - 1] * L.temp[j]
    result[j] = (rhs[j] - L.lowerDiagonal[j - 1] * result[j - 1]) / bet
  end

  for j = L.n-1:-1:2
    result[j] -= L.temp[j + 1] * result[j + 1]
  end

  result[1] -= L.temp[2] * result[2]

  return L, result
end
