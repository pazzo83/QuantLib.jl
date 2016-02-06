# Observer pattern (can be an observer or something that is observable)
type ObserverMixin
  observers::Vector
  observables::Vector
end

ObserverMixin() = ObserverMixin([], [])

# observable patterns
function add_observer!{T}(ob::Observer, obsv::T)
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
