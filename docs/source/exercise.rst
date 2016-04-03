Exercise
========

Types for option exercises

Two abstract types govern these structures:

.. code-block:: julia

    abstract Exercise
    abstract EarlyExercise <: Exercise


American Exercise
-----------------

An American option can be exercised at any time between two predefined dates; the first date might be omitted, in which case the option can be exercised at any time before the expiry.

.. code-block:: julia

    type AmericanExercise <: EarlyExercise
      dates::Vector{Date}
    end


.. function:: AmericanExercise(d1::Date, d2::Date)

    Constructor for an American Exercise with earliest date and latest date for exercise


Bermudan Exercise
-----------------

A Bermudan option can only be exercised at a set of fixed dates.

.. code-block:: juila

    type BermudanExercise <: EarlyExercise
      dates::Vector{Date}
    end


European Exercise
-----------------

A European option can only be exercised at one (expiry) date.

.. code-block:: julia

    type EuropeanExercise <: Exercise
      dates::Vector{Date}
    end

.. function:: EuropeanExercise(d::Date) = EuropeanExercise([d])

    Constructs a European Exercise with one expiration date


Rebated Exercise
---------------

In case of exercise the holder receives a rebate (if positive) or pays it (if negative) on the rebate settlement date.

.. code-block:: julia

    type RebatedExercise <: Exercise
      exercise::Exercise
      rebate::Float64
      rebateSettlementDays::Int
      rebatePaymentCalendar::BusinessCalendar
      rebatePaymentConvention::BusinessDayConvention
    end

.. function:: RebatedExercise(exercise::Exercise, rebate::Float64 = 0.0, rebateSettlementDays::Int = 0, rebatePaymentCalendar::BusinessCalendar = NullCalendar(), rebatePaymentConvention::BusinessDayConvention = Following())

    Default constructor for a RebatedExercise
