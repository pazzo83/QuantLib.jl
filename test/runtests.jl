using QuantLib

tests = ["cash_flows",
        "indexes",
        "bonds",
        "forwards",
        "options",
        "swaps",
        "swaptions",
        "fd_mesher"]

println("Running tests:")

for t in tests
    tfile = string(t, ".jl")
    println(" * $(tfile) ...")
    include(tfile)
end
