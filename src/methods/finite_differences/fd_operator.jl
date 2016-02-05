type TripleBandLinearOp{I <: Integer, FM <: FdmMesher}
  direction::I
  mesher::FM
  i0::Vector{I}
  i2::Vector{I}
  reverseIndex::Vector{I}
  lower::Vector{Float64}
  _diag::Vector{Float64}
  upper::Vector{Float64}
end

function TripleBandLinearOp(direction::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i0 = Vector{Int}(sz)
  i2 = Vector{Int}(sz)
  reverseIndex = Vector{Int}(sz)
  lower = Vector{Float64}(sz)
  _diag = Vector{Float64}(sz)
  upper = Vector{Float64}(sz)

  newDim = copy(mesher.layout.dim)
  newDim[1], newDim[direction] = newDim[direction], newDim[1] # swap

  newSpacing = FdmLinearOpLayout(newDim).spacing
  newSpacing[1], newSpacing[direction] = newSpacing[direction], newSpacing[1] # swap

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    ## TripleBandLinearOp part
    i0[i] = neighborhood(mesher.layout, i, coords, direction, -1)
    i2[i] = neighborhood(mesher.layout, i, coords, direction, 1)

    newIndex = dot(coords - 1, newSpacing) + 1

    reverseIndex[newIndex] = i

    iter_coords!(coords, mesher.layout.dim)
  end

  return TripleBandLinearOp(direction, mesher, i0, i2, reverseIndex, lower, _diag, upper)
end

function FirstDerivativeOp(direction::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i0 = Vector{Int}(sz)
  i2 = Vector{Int}(sz)
  reverseIndex = Vector{Int}(sz)
  lower = Vector{Float64}(sz)
  _diag = Vector{Float64}(sz)
  upper = Vector{Float64}(sz)

  newDim = copy(mesher.layout.dim)
  newDim[1], newDim[direction] = newDim[direction], newDim[1] # swap

  newSpacing = FdmLinearOpLayout(newDim).spacing
  newSpacing[1], newSpacing[direction] = newSpacing[direction], newSpacing[1] # swap

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    ## TripleBandLinearOp part
    i0[i] = neighborhood(mesher.layout, i, coords, direction, -1)
    i2[i] = neighborhood(mesher.layout, i, coords, direction, 1)

    newIndex = dot(coords - 1, newSpacing) + 1

    reverseIndex[newIndex] = i

    ## FirstDerivativeOp part
    hm = get_dminus(mesher, coords, direction)
    hp = get_dplus(mesher, coords, direction)

    zetam1 = hm * (hm + hp)
    zeta0 = hm * hp
    zetap1 = hp * (hm + hp)

    if coords[direction] == 1
      # upwinding scheme
      lower[i] = 0.0
      upper[i] = 1.0 / hp
      _diag[i] = -(upper[i])
    elseif coords[direction] == mesher.layout.dim[direction]
      # downwinding scheme
      _diag[i] = 1.0 / hm
      lower[i] = -_diag[i]
      upper[i] = 0.0
    else
      lower[i] = -hp / zetam1
      _diag[i] = (hp - hm) / zeta0
      upper[i] = hm / zetap1
    end
    iter_coords!(coords, mesher.layout.dim)
  end

  return TripleBandLinearOp(direction, mesher, i0, i2, reverseIndex, lower, _diag, upper)
end

function SecondDerivativeOp(direction::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i0 = Vector{Int}(sz)
  i2 = Vector{Int}(sz)
  reverseIndex = Vector{Int}(sz)
  lower = Vector{Float64}(sz)
  _diag = Vector{Float64}(sz)
  upper = Vector{Float64}(sz)

  newDim = copy(mesher.layout.dim)
  newDim[1], newDim[direction] = newDim[direction], newDim[1] # swap

  newSpacing = FdmLinearOpLayout(newDim).spacing
  newSpacing[1], newSpacing[direction] = newSpacing[direction], newSpacing[1] # swap

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    ## TripleBandLinearOp part
    i0[i] = neighborhood(mesher.layout, i, coords, direction, -1)
    i2[i] = neighborhood(mesher.layout, i, coords, direction, 1)

    newIndex = dot(coords - 1, newSpacing) + 1

    reverseIndex[newIndex] = i

    ## SecondDerivativeOp part
    hm = get_dminus(mesher, coords, direction)
    hp = get_dplus(mesher, coords, direction)

    zetam1 = hm * (hm + hp)
    zeta0 = hm * hp
    zetap1 = hp * (hm + hp)

    co = coords[direction]

    if co == 1 || co == mesher.layout.dim[direction]
      lower[i] = upper[i] = _diag[i] = 0.0
    else
      lower[i] = 2.0 / zetam1
      _diag[i] = -2.0 / zeta0
      upper[i] = 2.0 / zetap1
    end
    iter_coords!(coords, mesher.layout.dim)
  end

  return TripleBandLinearOp(direction, mesher, i0, i2, reverseIndex, lower, _diag, upper)
end

function mult!{T <: Number}(trpBandLinOp::TripleBandLinearOp, u::Vector{T})
  # probably should protect for array dims here
  trpBandLinOp.lower = trpBandLinOp.lower .* u
  trpBandLinOp._diag = trpBandLinOp._diag .* u
  trpBandLinOp.upper = trpBandLinOp.upper .* u

  return trpBandLinOp
end

function add!(trpBandLinOp::TripleBandLinearOp, m::TripleBandLinearOp)
  trpBandLinOp.lower += m.lower
  trpBandLinOp._diag += m._diag
  trpBandLinOp.upper += m.upper

  return trpBandLinOp
end

function axpyb!(trpBandLinOp::TripleBandLinearOp, a::Vector{Float64}, x::TripleBandLinearOp, y::TripleBandLinearOp, b::Vector{Float64})
  sz = trpBandLinOp.mesher.layout.size

  if length(a) == 0
    if length(b) == 0
      trpBandLinOp._diag[1:sz] = y._diag
      trpBandLinOp.lower[1:sz] = y.lower
      trpBandLinOp.upper[1:sz] = y.upper
    else
      addB = length(b) > 0 ? b[1:sz] : b[1]
      trpBandLinOp._diag[1:sz] = y._diag + addB
      trpBandLinOp.lower[1:sz] = y.lower
      trpBandLinOp.upper[1:sz] = y.upper
    end
  elseif length(b) == 0
    # this might be improved
    trpBandLinOp._diag[1:sz] = y._diag + (length(a) > 0 ? a .* x._diag : a[1] * x._diag)
    trpBandLinOp.lower[1:sz] = y.lower + (length(a) > 0 ? a .* x.lower : a[1] * x.lower)
    trpBandLinOp.upper[1:sz] = y.upper + (length(a) > 0 ? a .* x.upper : a[1] * x.upper)
  else
    addB = length(b) > 0 ? b[1:sz] : b[1]
    trpBandLinOp._diag[1:sz] = y._diag + (length(a) > 0 ? a .* x._diag : a[1] * x._diag) + addB
    trpBandLinOp.lower[1:sz] = y.lower + (length(a) > 0 ? a .* x.lower : a[1] * x.lower)
    trpBandLinOp.upper[1:sz] = y.upper + (length(a) > 0 ? a .* x.upper : a[1] * x.upper)
  end

  return trpBandLinOp
end

function apply(trpBandLinOp::TripleBandLinearOp, r::Vector{Float64})
  idx = trpBandLinOp.mesher.layout
  length(r) == idx.size || error("inconsistent length of r")

  retArray = Vector{Float64}(length(r))

  for i = 1:idx.size
    retArray[i] = r[trpBandLinOp.i0[i]] * trpBandLinOp.lower[i] + r[i] * trpBandLinOp._diag[i] + r[trpBandLinOp.i2[i]] * trpBandLinOp.upper[i]
  end

  return retArray
end

function solve_splitting(trpBandLinOp::TripleBandLinearOp, r::Vector{Float64}, a::Float64, b::Float64)
  layout = trpBandLinOp.mesher.layout
  length(r) == layout.size || error("inconsistent length of r")

  retArray = Vector{Float64}(length(r))
  tmp = Vector{Float64}(length(r))

  # solving a tridiagonal system (we could use Julia built in f'ns here)
  rim1 = trpBandLinOp.reverseIndex[1]
  bet = 1.0 / (a * trpBandLinOp._diag[rim1] + b)
  bet != 0.0 || error("division by zero")

  retArray[rim1] = r[rim1] * bet

  for j = 2:layout.size
    ri = trpBandLinOp.reverseIndex[j]
    tmp[j] = a * trpBandLinOp.upper[rim1] * bet

    bet = b + a * (trpBandLinOp._diag[ri] - tmp[j] * trpBandLinOp.lower[ri])
    bet != 0.0 || error("division by zero")

    bet = 1.0 / bet

    retArray[ri] = (r[ri] - a * trpBandLinOp.lower[ri] * retArray[rim1]) * bet
    rim1 = ri
  end

  for j = layout.size - 1:-1:2
    retArray[trpBandLinOp.reverseIndex[j]] -= tmp[j + 1] * retArray[trpBandLinOp.reverseIndex[j + 1]]
  end

  retArray[trpBandLinOp.reverseIndex[1]] -= tmp[2] * retArray[trpBandLinOp.reverseIndex[2]]

  return retArray
end

type SecondOrderMixedDerivativeOp{I <: Integer, FD <: FdmMesher} <: NinePointLinearOp
  d1::I
  d2::I
  i00::Vector{I}
  i10::Vector{I}
  i20::Vector{I}
  i01::Vector{I}
  i21::Vector{I}
  i02::Vector{I}
  i12::Vector{I}
  i22::Vector{I}
  a00::Vector{Float64}
  a10::Vector{Float64}
  a20::Vector{Float64}
  a01::Vector{Float64}
  a11::Vector{Float64}
  a21::Vector{Float64}
  a02::Vector{Float64}
  a12::Vector{Float64}
  a22::Vector{Float64}
  mesher::FD
end

function SecondOrderMixedDerivativeOp(d1::Int, d2::Int, mesher::FdmMesher)
  sz = mesher.layout.size
  i00 = Vector{Int}(sz)
  i10 = Vector{Int}(sz)
  i20 = Vector{Int}(sz)
  i01 = Vector{Int}(sz)
  i21 = Vector{Int}(sz)
  i02 = Vector{Int}(sz)
  i12 = Vector{Int}(sz)
  i22 = Vector{Int}(sz)
  a00 = Vector{Float64}(sz)
  a10 = Vector{Float64}(sz)
  a20 = Vector{Float64}(sz)
  a01 = Vector{Float64}(sz)
  a11 = Vector{Float64}(sz)
  a21 = Vector{Float64}(sz)
  a02 = Vector{Float64}(sz)
  a12 = Vector{Float64}(sz)
  a22 = Vector{Float64}(sz)

  coords = ones(Int, length(mesher.layout.dim))
  for i = 1:sz
    # NinePointLinearOp part
    i10[i] = neighborhood(mesher.layout, i, coords, d2, -1)
    i01[i] = neighborhood(mesher.layout, i, coords, d1, -1)
    i21[i] = neighborhood(mesher.layout, i, coords, d1,  1)
    i12[i] = neighborhood(mesher.layout, i, coords, d2,  1)
    i00[i] = neighborhood(mesher.layout, i, coords, d1, -1, d2, -1)
    i20[i] = neighborhood(mesher.layout, i, coords, d1,  1, d2, -1)
    i02[i] = neighborhood(mesher.layout, i, coords, d1, -1, d2,  1)
    i22[i] = neighborhood(mesher.layout, i, coords, d1,  1, d2,  1)

    # SecondOrderMixedDerivativeOp part
    hm_d1 = get_dminus(mesher, coords, d1)
    hp_d1 = get_dplus(mesher, coords, d1)
    hm_d2 = get_dminus(mesher, coords, d2)
    hp_d2 = get_dplus(mesher, coords, d2)

    zetam1 = hm_d1 * (hm_d1 + hp_d1)
    zeta0 = hm_d1 * hp_d1
    zetap1 = hp_d1 * (hm_d1 + hp_d1)

    phim1 = hm_d2 * (hm_d2 + hp_d2)
    phi0 = hm_d2 * hp_d2
    phip1 = hp_d2 * (hm_d2 + hp_d2)

    c1 = coords[d1]
    c2 = coords[d2]

    if c1 == 1 && c2 == 1
      # lower left corner
      a00[i] = a01[i] = a02[i] = a10[i] = a20[i] = 0.0
      a11[i] = a22[i] = 1.0 / (hp_d1 * hp_d2)
      a21[i] = a12[i] = -a11[i]
    elseif c1 == mesher.layout.dim[d1] && c2 == 1
      # upper left corner
      a22[i] = a21[i] = a20[i] = a10[i] = a00[i] = 0.0
      a01[i] = a12[i] = 1.0 / (hm_d1 * hp_d2)
      a11[i] = a02[i] = -a01[i]
    elseif c1 == 1 && c2 == mesher.layout.dim[d2]
      # lower right corner
      a00[i] = a01[i] = a02[i] = a12[i] = a22[i] = 0.0
      a10[i] = a21[i] = 1.0 / (hp_d1 * hm_d2)
      a20[i] = a11[i] = -a10[i]
    elseif c1 == mesher.layout.dim[d1] && c2 == mesher.layout.dim[d2]
      # upper right corner
      a20[i] = a21[i] = a22[i] = a12[i] = a02[i] = 0.0
      a00[i] = a11[i] = 1.0 / (hm_d1 * hm_d2)
      a10[i] = a01[i] = -a00[i]
    elseif c1 == 1
      # lower side
      a00[i] = a01[i] = a02[i] = 0.0
      a10[i] = hp_d2 / (hp_d1 * phim1)
      a20[i] = -a10[i]
      a21[i] = (hp_d2 - hm_d2) / (hp_d1 * phi0)
      a11[i] = -a21[i]
      a22[i] = hm_d2 / (hp_d1 * phip1)
      a12[i] = -a22[i]
    elseif c1 == mesher.layout.dim[d1]
      # upper side
      a20[i] = a21[i] = a22[i] = 0.0
      a00[i] = hp_d2 / (hm_d1 * phim1)
      a10[i] = -a00[i]
      a11[i] = (hp_d2 - hm_d2) / (hm_d1 * phi0)
      a01[i] = -a11[i]
      a12[i] = hm_d2 / (hm_d1 * phip1)
      a02[i] = -a12[i]
    elseif c2 == 1
      # left side
      a00[i] = a10[i] = a20[i] = 0.0
      a01[i] = hp_d1 / (zetam1 * hp_d2)
      a02[i] = -a01[i]
      a12[i] = (hp_d1 - hm_d1) / (zeta0 * hp_d2)
      a11[i] = -a12[i]
      a22[i] = hm_d1 / (zetap1 * hp_d2)
      a21[i] = -a22[i]
    elseif c2 == mesher.layout.dim[d2]
      # right side
      a22[i] = a12[i] = a02[i] = 0.0
      a00[i] = hp_d1 / (zetam1 * hm_d2)
      a01[i] = -a00[i]
      a11[i] = (hp_d1 - hm_d1) / (zeta0 * hm_d2)
      a10[i] = -a11[i]
      a21[i] = hm_d1 / (zetap1 * hm_d2)
      a20[i] = -a21[i]
    else
      a00[i] = hp_d1 * hp_d2 / (zetam1 * phim1)
      a10[i] = -(hp_d1 - hm_d1) * hp_d2 / (zeta0 * phim1)
      a20[i] = -hm_d1 * hp_d2 / (zetap1 * phim1)
      a01[i] = -hp_d1 * (hp_d2 - hm_d2) / (zetam1 * phi0)
      a11[i] = (hp_d1 - hm_d1) * (hp_d2 - hm_d2) / (zeta0 * phi0)
      a21[i] = hm_d1 * (hp_d2 - hm_d2) / (zetap1 * phi0)
      a02[i] = -hp_d1 * hm_d2 / (zetam1 * phip1)
      a12[i] = hm_d2 * (hp_d1 - hm_d1) / (zeta0 * phip1)
      a22[i] = hm_d1 * hm_d2 / (zetap1 * phip1)
    end
    iter_coords!(coords, mesher.layout.dim)
  end

  return SecondOrderMixedDerivativeOp(d1, d2, i00, i10, i20, i01, i21, i02, i12, i22, a00, a10, a20, a01, a11, a21, a02, a12, a22, mesher)
end

function mult!{T <: Number}(ninePointLin::NinePointLinearOp, u::Vector{T})
  ninePointLin.a11 = ninePointLin.a11 .* u
  ninePointLin.a01 = ninePointLin.a01 .* u
  ninePointLin.a10 = ninePointLin.a10 .* u
  ninePointLin.a21 = ninePointLin.a21 .* u
  ninePointLin.a22 = ninePointLin.a22 .* u
  ninePointLin.a00 = ninePointLin.a00 .* u
  ninePointLin.a02 = ninePointLin.a02 .* u
  ninePointLin.a20 = ninePointLin.a20 .* u
  ninePointLin.a12 = ninePointLin.a12 .* u

  return ninePointLin
end

function apply(ninePointLin::NinePointLinearOp, u::Vector{Float64})
  length(u) == ninePointLin.mesher.layout.size || error("inconsistent length of r")

  retVal = ninePointLin.a00 .* u[ninePointLin.i00] + ninePointLin.a01 .* u[ninePointLin.i01] + ninePointLin.a02 .* u[ninePointLin.i02] + ninePointLin.a10 .* u[ninePointLin.i10] +
          ninePointLin.a11 .* u + ninePointLin.a12 .* u[ninePointLin.i12] + ninePointLin.a20 .* u[ninePointLin.i20] + ninePointLin.a21 .* u[ninePointLin.i21] +
          ninePointLin.a22 .* u[ninePointLin.i22]

  return retVal
end

type FdmG2Op{I <: Integer} <: FdmLinearOpComposite
  direction1::I
  direction2::I
  x::Vector{Float64}
  y::Vector{Float64}
  dxMap::TripleBandLinearOp
  dyMap::TripleBandLinearOp
  corrMap::SecondOrderMixedDerivativeOp
  mapX::TripleBandLinearOp
  mapY::TripleBandLinearOp
  model::G2
end

function FdmG2Op(mesher::FdmMesher, model::G2, direction1::Int, direction2::Int)
  x = get_locations(mesher, direction1)
  y = get_locations(mesher, direction2)

  dxMap = add!(mult!(FirstDerivativeOp(direction1, mesher), (-x * get_a(model))),
              mult!(SecondDerivativeOp(direction1, mesher), (0.5 * get_sigma(model) * get_sigma(model) * ones(mesher.layout.size))))

  dyMap = add!(mult!(FirstDerivativeOp(direction2, mesher), (-y * get_b(model))),
              mult!(SecondDerivativeOp(direction2, mesher), (0.5 * get_eta(model) * get_eta(model) * ones(mesher.layout.size))))

  corrMap = mult!(SecondOrderMixedDerivativeOp(direction1, direction2, mesher), fill(get_rho(model) * get_sigma(model) * get_eta(model), mesher.layout.size))

  mapX = TripleBandLinearOp(direction1, mesher)
  mapY = TripleBandLinearOp(direction2, mesher)

  return FdmG2Op(direction1, direction2, x, y, dxMap, dyMap, corrMap, mapX, mapY, model)
end

type FdmHullWhiteOp{I <: Integer} <: FdmLinearOpComposite
  direction::I
  x::Vector{Float64}
  dzMap::TripleBandLinearOp
  mapT::TripleBandLinearOp
  model::HullWhite
end

function FdmHullWhiteOp(mesher::FdmMesher, model::HullWhite, direction::Int)
  x = get_locations(mesher, direction)

  dzMap = add!(mult!(FirstDerivativeOp(direction, mesher), (-x * get_a(model))),
              mult!(SecondDerivativeOp(direction, mesher), (0.5 * get_sigma(model) * get_sigma(model)) * ones(mesher.layout.size)))

  mapT = TripleBandLinearOp(direction, mesher)

  return FdmHullWhiteOp(direction, x, dzMap, mapT, model)
end

# OP methods #
function set_time!(op::FdmG2Op, t1::Float64, t2::Float64)
  dynamics = get_dynamics(op.model)

  phi = 0.5 * (short_rate(dynamics, t1, 0.0, 0.0) + short_rate(dynamics, t2, 0.0, 0.0))

  hr = -0.5 * (op.x + op.y + phi)
  axpyb!(op.mapX, Vector{Float64}(0), op.dxMap, op.dxMap, hr)
  axpyb!(op.mapY, Vector{Float64}(0), op.dyMap, op.dyMap, hr)

  return op
end

function set_time!(op::FdmHullWhiteOp, t1::Float64, t2::Float64)
  dynamics = get_dynamics(op.model)

  phi = 0.5 * (short_rate(dynamics, t1, 0.0) + short_rate(dynamics, t2, 0.0))

  axpyb!(op.mapT, Vector{Float64}(0), op.dzMap, op.dzMap, -(op.x + phi))

  return op
end

function apply_direction(op::FdmG2Op, direction::Int, r::Vector{Float64})
  if direction == op.direction1
    return apply(op.mapX, r)
  elseif direction == op.direction2
    return apply(op.mapY, r)
  else
    return zeros(length(r))
  end
end

function apply_direction(op::FdmHullWhiteOp, direction::Int, r::Vector{Float64})
  if direction == op.direction
    return apply(op.mapT, r)
  else
    return zeros(length(r))
  end
end

function solve_splitting(op::FdmG2Op, direction::Int, r::Vector{Float64}, a::Float64)
  if direction == op.direction1
    return solve_splitting(op.mapX, r, a, 1.0)
  elseif direction == op.direction2
    return solve_splitting(op.mapY, r, a, 1.0)
  else
    return zeros(length(r))
  end
end

function solve_splitting(op::FdmHullWhiteOp, direction::Int, r::Vector{Float64}, a::Float64)
  if direction == op.direction
    return solve_splitting(op.mapT, r, a, 1.0)
  else
    return zeros(length(r))
  end
end

apply_mixed(op::FdmG2Op, r::Vector{Float64}) = apply(op.corrMap, r)
apply_mixed(op::FdmHullWhiteOp, r::Vector{Float64}) = zeros(length(r))

apply(op::FdmG2Op, r::Vector{Float64}) = apply(op.mapX, r) + apply(op.mapY, r) + apply_mixed(op, r)
apply(op::FdmHullWhiteOp, r::Vector{Float64}) = apply(op.mapT, r)

get_size(op::FdmG2Op) = 2
get_size(op::FdmHullWhiteOp) = 1
