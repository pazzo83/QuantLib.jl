# Basic model stuff
# model types
type AffineModelType <: ModelType end
type TermStructureConsistentModelType <: ModelType end

# type aliases
typealias AffineModel Model{AffineModelType}
typealias TermStructureConsistentModel Model{TermStructureConsistentModelType}

type ModelCommon
  observe::ObserverMixin
end

ModelCommon() = ModelCommon(ObserverMixin())

get_observers(m::Model) = m.common.observe.observers
