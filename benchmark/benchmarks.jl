#!/usr/bin/env julia
# benchmark/benchmarks.jl — BenchmarkTools suite for SNA.jl's hot loops.
#
# Covers the edge-driven (Batagelj–Mrvar) triad census that replaced the
# O(n³) triple loop: on sparse networks its cost is O(m · d̄), so with the
# mean degree held fixed the census scales ~linearly in n rather than
# cubically. Benchmarked at n = 500 and n = 2000, directed and undirected.
#
# Defines the standard `SUITE::BenchmarkGroup`. Run standalone with
#     julia --project=benchmark benchmark/benchmarks.jl
# which tunes + runs the suite and prints one tab-separated `BENCHJL` line
# per benchmark (consumed by the site repo's tools/run_benchmarks.jl).

using BenchmarkTools
using Network
using Random
using SNA

# ---------------------------------------------------------------------------
# Fixtures: sparse Erdős–Rényi networks with fixed expected mean degree
# ---------------------------------------------------------------------------

const MEAN_DEGREE = 10
const SIZES = (500, 2000)

function er_network(rng::AbstractRNG, n::Int; directed::Bool=false)
    net = network(n; directed=directed)
    m = directed ? MEAN_DEGREE * n : (MEAN_DEGREE * n) ÷ 2
    while ne(net) < m
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i == j && continue
        add_edge!(net, i, j)
    end
    return net
end

const NETS_DIR = Dict(n => er_network(Random.Xoshiro(n), n; directed=true)
                      for n in SIZES)
const NETS_UND = Dict(n => er_network(Random.Xoshiro(n + 1), n; directed=false)
                      for n in SIZES)

# ---------------------------------------------------------------------------
# Suite
# ---------------------------------------------------------------------------

const SUITE = BenchmarkGroup()

let g = addgroup!(SUITE, "triad_census")
    for n in SIZES
        g["directed_n$(n)"] = @benchmarkable triad_census($(NETS_DIR[n]))
        g["undirected_n$(n)"] = @benchmarkable triad_census($(NETS_UND[n]))
    end
end

# ---------------------------------------------------------------------------
# Standalone entry point
# ---------------------------------------------------------------------------

function print_benchjl(results::BenchmarkGroup)
    for (path, trial) in BenchmarkTools.leaves(results)
        est = median(trial)
        println("BENCHJL\t", join(path, "/"), "\t",
                BenchmarkTools.time(est), "\t",
                BenchmarkTools.allocs(est), "\t",
                BenchmarkTools.memory(est))
    end
end

function main()
    tune!(SUITE)
    results = run(SUITE; verbose=false, seconds=1)
    print_benchjl(results)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
