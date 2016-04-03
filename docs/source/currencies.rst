Currencies
==========

Derived from Ito.jl, the Currency type contains a variety of currencies worldwide, which are generated at compile time.

.. code-block:: julia

    immutable Currency <: AbstractCurrency
      name::AbstractString
      code::AbstractString
      numeric::Int
      symbol::AbstractString
      fractionSymbol::AbstractString
      fractionsPerUnit::Int
      rounding::Function
      formatString::AbstractString
    end

Example of usage:

.. code-block:: julia

    EURCurrency() # Euro
    USDCurrency() # USD
