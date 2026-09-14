using BenchmarkTools
using Networks
using Random
using SNA
using Test

# Fixed mean degree isolates the census scaling from graph density. Build
# fixtures outside the measured call, with independent seeded RNGs.
function sparse_fixture(n)
    rng = Xoshiro(n)
    net = Network(n)
    while ne(net) < 6n
        i, j = rand(rng, 1:n), rand(rng, 1:n)
        i != j && add_edge!(net, i, j)
    end
    return net
end

@testset "Sparse triad census scaling" begin
    small, large = sparse_fixture(400), sparse_fixture(1600)
    triad_census(small); triad_census(large)
    a_small = @allocated triad_census(small)
    a_large = @allocated triad_census(large)
    @test a_large < 6a_small
    t_small = @belapsed triad_census($small) seconds=0.3
    t_large = @belapsed triad_census($large) seconds=0.3
    # Four times as many actors/edges should remain near linear, with enough
    # headroom for CI timing noise; an O(n³) census would grow by 64×.
    @test t_large < 10t_small
    @test sum(triad_census(large)) == binomial(nv(large), 3)
    println("SCALING\ttriad_census\ttime_ratio=", t_large / t_small,
            "\tallocation_ratio=", a_large / a_small)
end
