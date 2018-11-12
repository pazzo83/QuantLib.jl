using Test
using QuantLib

mesher1 = Concentrating1dMesher(0.1, 2.5, 1001)
@test mesher1.locations[1] == range(0.1, stop=2.5, length=1001)[1]
mesher2 = Concentrating1dMesher(0.1, 2.5, 1001, (1.1, 0.25), true)
@test mesher2.locations[2] ≈ 0.10333775170976389 atol=1.0e-8