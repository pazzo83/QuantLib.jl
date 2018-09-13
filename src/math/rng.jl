using StatsFuns
using Sobol
using Random

abstract type AbstractRandomSequenceGenerator end
mutable struct PseudoRandomRSG <: AbstractRandomSequenceGenerator
  rng::MersenneTwister
  dimension::Int
  values::Vector{Float64}
  weight::Float64
end

PseudoRandomRSG(seed::Int, dimension::Int = 1, weight::Float64 = 1.0) = PseudoRandomRSG(MersenneTwister(seed), dimension, zeros(dimension), weight)

mutable struct InverseCumulativeRSG <: AbstractRandomSequenceGenerator
  rng::MersenneTwister
  dimension::Int
  values::Vector{Float64}
  weight::Float64
end

InverseCumulativeRSG(seed::Int, dimension::Int = 1, weight::Float64 = 1.0) = InverseCumulativeRSG(MersenneTwister(seed), dimension, zeros(dimension), weight)

mutable struct SobolRSG <: AbstractRandomSequenceGenerator
  rng::Sobol.SobolSeq
  dimension::Int
  values::Vector{Float64}
  weight::Float64
end

SobolRSG(dimension::Int = 1, weight::Float64 = 1.0) = SobolRSG(SobolSeq(dimension), dimension, zeros(dimension), weight)

mutable struct SobolInverseCumulativeRSG <: AbstractRandomSequenceGenerator
  rng::Sobol.SobolSeq
  dimension::Int
  values::Vector{Float64}
  weight::Float64
end

SobolInverseCumulativeRSG(dimension::Int = 1, weight::Float64 = 1.0) = SobolInverseCumulativeRSG(SobolSeq(dimension), dimension, zeros(dimension), weight)

function next_sequence!(rsg::InverseCumulativeRSG)
  # we can probably use map here with the norminvcdf
  # like this: map!(norminvcdf, rsg.values, rand(rsg, length(rsg.values)))
  @inbounds @simd for i in eachindex(rsg.values)
    x = rand(rsg.rng) # get random number
    rsg.values[i] = norminvcdf(x)
  end

  return rsg.values, rsg.weight
end

function next_sequence!(rsg::SobolRSG)
  # we can probably use map here with the norminvcdf
  # like this: map!(norminvcdf, rsg.values, rand(rsg, length(rsg.values)))
  rsg.values = next(rsg.rng)
  # for i in eachindex(rsg.values)
  #   x = next(rsg.rng) # get random number
  #   rsg.values[i] = x
  # end

  return rsg.values, rsg.weight
end

function next_sequence!(rsg::SobolInverseCumulativeRSG)
  seq = next(rsg.rng)

  @simd for i in eachindex(rsg.values)
    @inbounds rsg.values[i] = norminvcdf(seq[i])
  end

  return rsg.values, rsg.weight
end

function next_sequence!(rsg::PseudoRandomRSG)
  rsg.values = rand(rsg.rng, rsg.dimension)
  return rsg.values, rsg.weight
end

last_sequence(rsg::AbstractRandomSequenceGenerator) = rsg.values, rsg.weight

function init_sequence_generator!(rsg::SobolRSG, dimension::Int)
  rsg.dimension = dimension
  rsg.values = zeros(dimension)
  rsg.rng = SobolSeq(dimension)

  return rsg
end

function init_sequence_generator!(rsg::AbstractRandomSequenceGenerator, dimension::Int)
  # MersenneTwister already init-ed
  rsg.dimension = dimension
  rsg.values = zeros(dimension)

  return rsg
end
