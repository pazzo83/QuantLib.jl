# Basic model stuff
# model types
type AffineModelType <: ModelType end
type TermStructureConsistentModelType <: ModelType end
type CalibratedModelType <: ModelType end # general model type

# type aliases
typealias AffineModel Model{AffineModelType}
typealias TermStructureConsistentModel Model{TermStructureConsistentModelType}
typealias CalibratedModel Model{CalibratedModelType}

type ModelCommon
  observe::ObserverMixin
end

ModelCommon() = ModelCommon(ObserverMixin())

get_observers(m::Model) = m.common.observe.observers
