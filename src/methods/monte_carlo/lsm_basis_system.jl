type Monomial <: LsmBasisSystemPolynomType end

type MonomialFunction{I <: Integer} <: LSMBasisSystemFunction
  order::I
end

function call(m::MonomialFunction, x::Float64)
  ret = 1.0
  for i = 1:m.order
    ret *= x
  end

  return ret
end

get_type(::Monomial) = MonomialFunction{Int}

function path_basis_system!(::Monomial, order::Int, v::Vector)
  for i = 1:order + 1
    v[i] = MonomialFunction(i - 1)
  end

  return v
end
