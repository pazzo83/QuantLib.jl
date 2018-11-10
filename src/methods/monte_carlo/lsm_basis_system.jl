struct Monomial <: LsmBasisSystemPolynomType end

struct MonomialFunction <: LSMBasisSystemFunction
  order::Int
end

function (m::MonomialFunction)(x::Float64)
  ret = 1.0
  @simd for i = 1:m.order
    ret *= x
  end

  return ret
end

get_type(::Monomial) = MonomialFunction{Int}

function path_basis_system!(::Monomial, order::Int, v::Vector)
  @simd for i = 1:order + 1
    @inbounds v[i] = MonomialFunction(i - 1)
  end

  return v
end
