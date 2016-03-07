using StatsFuns

abstract AbstractRandomSequenceGenerator
type InverseCumulativeRSG{I <: Integer} <: AbstractRandomSequenceGenerator
  rng::MersenneTwister
  dimension::I
  values::Vector{Float64}
  weight::Float64
end

InverseCumulativeRSG(seed::Int, dimension::Int = 1, weight::Float64 = 1.0) = InverseCumulativeRSG(MersenneTwister(seed), dimension, zeros(dimension), weight)

function next_sequence!(rsg::InverseCumulativeRSG)
  # we can probably use map here with the norminvcdf
  # like this: map!(norminvcdf, rsg.values, rand(rsg, length(rsg.values)))
  for i in eachindex(rsg.values)
    x = rand(rsg) # get random number
    rsg.values[i] = norminvcdf(x)
  end

  return rsg.values, rsg.weight
end

last_sequence(rsg::InverseCumulativeRSG) = rsg.values, rsg.weight
