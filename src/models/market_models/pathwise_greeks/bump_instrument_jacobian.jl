type VolatilityBumpInstrumentJacobianSwaption
  startIndex::Int
  endIndex::Int
end

type VolatilityBumpInstrumentJacobianCap
  startIndex::Int
  endIndex::Int
  strike::Float64
end

type VolatilityBumpInstrumentJacobian
  bumps::VegaBumpCollection
  swaptions::Vector{VolatilityBumpInstrumentJacobianSwaption}
  caps::Vector{VolatilityBumpInstrumentJacobianCap}
  computed::BitArray{1}
  allComputed::Bool
  derivatives::Vector{Vector{Float64}}
  onePercentBumps::Vector{Vector{Float64}}
  bumpMatrix::Matrix{Float64}
end

function VolatilityBumpInstrumentJacobian(bumps::VegaBumpCollection,
                                          swaptions::Vector{VolatilityBumpInstrumentJacobianSwaption},
                                          caps::Vector{VolatilityBumpInstrumentJacobianCap})

  swaptionsPlusCaps = length(swaptions) + length(caps)
  computed = falses(swaptionsPlusCaps)
  derivatives = [zeros(number_of_bumps(bumps)) for _ = 1:swaptionsPlusCaps]
  bumpMatrix = Matrix{Float64}(swaptionsPlusCaps, number_of_bumps(bumps))

  onePercentBumps = deepcopy(derivatives)
  allComputed = false

  return VolatilityBumpInstrumentJacobian(bumps, swaptions, caps, computed, allComputed, derivatives, onePercentBumps, bumpMatrix)
end

function derivatives_volatility!(voljacobian::VolatilityBumpInstrumentJacobian, j::Int)
  j <= length(voljacobian.swaptions) + length(voljacobian.caps) || error("too high index passed to derivatives_volatility")

  if voljacobian.computed[j]
    return voljacobian.derivatives[j]
  end

  resize!(voljacobian.derivatives[j], number_of_bumps(voljacobian.bumps))
  resize!(voljacobian.onePercentBumps[j], number_of_bumps(voljacobian.bumps))

  sizesq = 0.0
  voljacobian.computed[j] = true

  initj = j

  if j <= length(voljacobian.swaptions)
    # it's a swaption
    thisPseudo = SwaptionPseudoDerivative(associated_model(voljacobian.bumps), voljacobian.swaptions[j].startIndex, voljacobian.swaptions[j].endIndex)

    for k = 1:number_of_bumps(voljacobian.bumps)
      v = 0.0
      for i = voljacobian.allBumps[k].stepBegin:voljacobian.allBumps[k].stepEnd
        fullDerivative = thisPseudo.volatilityDerivatives[i]
        for f = voljacobian.allBumps[k].factorBegin:voljacobian.allBumps[k].factorEnd, r = voljacobian.allBumps[k].rateBegin:voljacobian.allBumps[k].rateEnd
          v += fullDerivative[r, f]
        end
      end

      voljacobian.derivatives[j][k] = v
      sizesq += v * v
    end
  else
    # it's a cap
    j -= length(voljacobian.swaptions) - 1 # need to get back to index 1
  end

  return voljacobian.derivatives[initj]
end

function get_all_one_percent_bumps!(voljacobian::VolatilityBumpInstrumentJacobian)
  if ~voljacobian.allComputed
    for i = 1:length(voljacobian.swaptions) + length(voljacobian.caps)
      derivatives_volatility!(voljacobian, i)
    end
  end

  voljacobian.allComputed = true

  return bumpMatrix
end


type OrthogonalizedBumpFinder
  derivativesProducer::VolatilityBumpInstrumentJacobian
  multiplierCutoff::Float64
  tolerance::Float64
end

OrthogonalizedBumpFinder(bumps::VegaBumpCollection,
                        swaptions::Vector{VolatilityBumpInstrumentJacobianSwaption},
                        caps::Vector{VolatilityBumpInstrumentJacobianCap},
                        multiplierCutoff::Float64,
                        tolerance::Float64) = OrthogonalizedBumpFinder(VolatilityBumpInstrumentJacobian(bumps, swaptions, caps), multiplierCutoff, tolerance)


function get_vega_bumps!(theBumps::Vector{Matrix{Float64}})
 # todo
end
