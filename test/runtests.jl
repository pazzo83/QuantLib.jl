using QuantLib

tests = ["cash_flows",
        "indexes",
        "bonds"]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end
