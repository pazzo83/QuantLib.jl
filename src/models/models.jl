# Basic model stuff
# model types
struct AffineModelType <: ModelType end
struct TermStructureConsistentModelType <: ModelType end
struct CalibratedModelType <: ModelType end # general model type

# type aliases
const AffineModel = Model{AffineModelType}
const TermStructureConsistentModel = Model{TermStructureConsistentModelType}
const CalibratedModel = Model{CalibratedModelType}

mutable struct ModelCommon
  observe::ObserverMixin
end

ModelCommon() = ModelCommon(ObserverMixin())

get_observers(m::Model) = m.common.observe.observers
