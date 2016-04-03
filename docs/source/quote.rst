Quotes
======

Basic Quote type for a price quote.

.. code-block:: julia

    type Quote
      value::Float64
      is_valid::Bool
    end

.. function:: Quote(value::Float64, is_valid::Bool = true) = new(value, is_valid)

    Default constructor for a Quote
