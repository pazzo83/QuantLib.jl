type SobolDiagonalOrdering <: SobolOrdering end
type SobolStepsOrdering <: SobolOrdering end
type SobolFactorsOrdering <: SobolOrdering end

type SobolBrownianGenerator{O <: SobolOrdering} <: BrownianGenerator
  factors::Int
  steps::Int
  ordering::O
  generator::SobolInverseCumulativeRSG
  bridge::BrownianBridge
  lastStep::Int
  orderedIndices::Vector{Vector{Int}}
  bridgedVariates::Vector{Vector{Float64}}
end

function SobolBrownianGenerator(factors::Int, steps::Int, ordering::SobolOrdering, ::Int)
  generator = SobolInverseCumulativeRSG(factors * steps)
  bridge = BrownianBridge(steps)
  orderedIndices = Vector{Int}[Vector{Int}(steps) for _ = 1:factors]
  bridgedVariates = Vector{Float64}[Vector{Float64}(steps) for _ = 1:factors]

  return SobolBrownianGenerator(factors, steps, ordering, generator, bridge, 0, orderedIndices, bridgedVariates)
end

type SobolBrownianGeneratorFactory{O <: SobolOrdering} <: BrownianGeneratorFactory
  ordering::O
  seed::Int
end

create(sob::SobolBrownianGeneratorFactory, factors::Int, steps::Int) = SobolBrownianGenerator(factors, steps, sob.ordering, sob.seed)
