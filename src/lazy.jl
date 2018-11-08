## Allows for lazy calculation (when requested)
## Also will have observer mixin in here

mutable struct LazyMixin
  calculated::Bool
  frozen::Bool
  observe::ObserverMixin

  LazyMixin() = new(false, false, ObserverMixin())
end

is_calculated(lazy::LazyObject) = lazy.lazyMixin.calculated
is_frozen(lazy::LazyObject) = lazy.lazyMixin.frozen

calculated!(lazy::LazyObject, setting::Bool = true) = lazy.lazyMixin.calculated = setting

function calculate!(lazy::LazyObject)
  if !is_calculated(lazy)
    calculated!(lazy)
    perform_calculations!(lazy)
  end

  return lazy
end

function recalculate!(lazy::LazyObject)
  calculated!(lazy, false)
  calculate!(lazy)

  return lazy
end

get_observers(lazy::LazyObject) = lazy.lazyMixin.observe.observers
