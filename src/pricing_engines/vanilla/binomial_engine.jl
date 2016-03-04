type BinomialVanillaEngine{P <: AbstractBlackScholesProcess, I <: Integer, T <: BinomialTreeType} <: AbstractVanillaEngine
  process::P
  timeSteps::I
  treeClass::T
end
