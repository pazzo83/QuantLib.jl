Welcome to QuantLib.jl's documentation!
=======================================
*QuantLib.jl* is a Julia package that provides a pure Julia version of the popular open-source quantitative finance library QuantLib.

This documentation is largely derived from QuantLib's documentation, with some alterations based on the Julia implementation.

**Contents:**

.. toctree::
   :maxdepth: 2

   getting_started.rst
   cashflows.rst
   currencies.rst
   indexes.rst
   instruments.rst
   interest_rates.rst
   math.rst
   methods.rst
   models.rst
   pricing_engines.rst
   process.rst
   quote.rst
   term_structures.rst
   time.rst


**QuantLib.jl Settings**

QuantLib.jl, upon loading, generates a Settings object which contains the global evaluation date.

.. code-block:: julia

    type Settings
      evaluation_date::Date
    end

.. function:: set_eval_date!(sett::Settings, d::Date)

    Update the global evaluation date



.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
