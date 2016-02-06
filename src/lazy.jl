## Allows for lazy calculation (when requested)
## Also will have observer mixin in here

type LazyMixin
  calculated::Bool
  frozen::Bool
  observe::ObserverMixin

  LazyMixin() = new(false, false, ObserverMixin())
end

is_calculated{L <: LazyObject}(lazy::L) = lazy.lazyMixin.calculated
is_frozen{L <: LazyObject}(lazy::L) = lazy.lazyMixin.frozen

calculated!{L <: LazyObject}(lazy::L, setting::Bool = true) = lazy.lazyMixin.calculated = setting

function calculate!{L <: LazyObject}(lazy::L)
  if !is_calculated(lazy)
    calculated!(lazy)
    perform_calculations!(lazy)
  end

  return lazy
end

function recalculate!{L <: LazyObject}(lazy::L)
  calculated!(lazy, false)
  calculate!(lazy)

  return lazy
end

get_observers(lazy::LazyObject) = lazy.lazyMixin.observe.observers
