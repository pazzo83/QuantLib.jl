type VegaBumpCluster
  factorBegin::Int
  factorEnd::Int
  rateBegin::Int
  rateEnd::Int
  stepBegin::Int
  stepEnd::Int
end

type VegaBumpCollection{M <: AbstractMarketModel}
  allBumps::Vector{VegaBumpCluster}
  associatedVolStructure::M
  checked::Bool
  nonOverlapped::Bool
  isFull::Bool
end

function VegaBumpCollection(volStructure::AbstractMarketModel, factorwiseBumping::Bool)
  steps = volStructure.numberOfSteps
  rates = volStructure.numberOfRates
  factors = volStructure.numberOfFactors

  allBumps = VegaBumpCluster[]

  for s = 1:steps
    for r = volStructure.evolution.firstAliveRate[s]:rates
      if factorwiseBumping
        for f = 1:factors
          push!(allBumps, VegaBumpCluster(f, f+1, r, r+1, s, s+1))
        end
      else
        push!(allBumps, VegaBumpCluster(1, factors+1, r, r+1, s, s+1))
      end
    end
  end

  return VegaBumpCollection(allBumps, volStructure, true, true, true)
end

number_of_bumps(vegaCol::VegaBumpCollection) = length(vegaCol.allBumps)
associated_model(vegaCol::VegaBumpCollection) = vegaCol.associatedVolStructure
