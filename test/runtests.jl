using QuantLib

tests = ["cash_flows",
        "indexes",
        "bonds",
        "forwards",
        "options",
        "swaps",
        "swaptions",
        "fd_mesher",
        "curves"]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end
