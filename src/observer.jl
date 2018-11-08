# Observer pattern (can be an observer or something that is observable)
mutable struct ObserverMixin
  observers::Vector
  observables::Vector
end

ObserverMixin() = ObserverMixin([], [])

# observable patterns
function add_observer!(ob::Observer, obsv::T) where {T}
  if ~in(obsv, get_observers(ob))
    push!(get_observers(ob), obsv)
  end

  return ob
end

function notify_observers!(ob::Observer)
  for obsv in get_observers(ob)
    update!(obsv)
  end

  return ob
end
