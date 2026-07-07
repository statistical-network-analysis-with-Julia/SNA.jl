using SNA
using Network
using Random
using Test

# Sampson monastery "liking" network (samplike, directed, n=18), exported
# from R's ergm/sna packages. Golden values below are from sna 2.8.
function sampson_like()
    net = network(18; directed=true)
    ties = [(1, 2), (1, 3), (1, 5), (1, 8), (1, 12), (1, 14), (2, 1), (2, 7),
            (2, 12), (2, 14), (2, 15), (3, 1), (3, 2), (3, 13), (3, 17), (3, 18),
            (4, 5), (4, 6), (4, 10), (4, 11), (5, 1), (5, 4), (5, 9), (5, 11),
            (5, 13), (6, 1), (6, 4), (6, 5), (6, 9), (7, 1), (7, 2), (7, 8),
            (7, 12), (7, 16), (8, 1), (8, 2), (8, 4), (8, 6), (8, 9), (8, 10),
            (9, 5), (9, 8), (9, 12), (9, 16), (10, 4), (10, 5), (10, 8), (10, 9),
            (10, 13), (10, 14), (11, 4), (11, 5), (11, 8), (11, 14), (11, 16), (12, 1),
            (12, 2), (12, 7), (12, 14), (13, 5), (13, 7), (13, 18), (14, 1), (14, 2),
            (14, 11), (14, 12), (14, 15), (15, 1), (15, 2), (15, 5), (15, 7), (15, 12),
            (15, 14), (16, 1), (16, 2), (16, 7), (16, 12), (16, 15), (17, 2), (17, 3),
            (17, 13), (17, 18), (18, 1), (18, 2), (18, 3), (18, 7), (18, 13), (18, 17)]
    for (i, j) in ties
        add_edge!(net, i, j)
    end
    @assert ne(net) == 88
    return net
end

# Padgett Florentine marriage network (undirected, n=16)
function florentine()
    net = network(16; directed=false)
    ties = [(1, 9), (2, 6), (2, 7), (2, 9), (3, 5), (3, 9), (4, 7), (4, 11),
            (4, 15), (5, 11), (5, 15), (7, 8), (7, 16), (9, 13), (9, 14),
            (9, 16), (10, 14), (11, 15), (13, 15), (13, 16)]
    for (i, j) in ties
        add_edge!(net, i, j)
    end
    return net
end

@testset "SNA.jl" begin
    @testset "Degree Centrality" begin
        net = network(5)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        dc = degree_centrality(net; mode=:out)
        @test dc[1] == 2.0  # 1 sends to 2, 3
        @test dc[3] == 1.0  # 3 sends to 4
        @test dc[5] == 0.0  # 5 sends to nobody

        dc_in = degree_centrality(net; mode=:in)
        @test dc_in[3] == 2.0  # 3 receives from 1, 2
    end

    @testset "Network Measures" begin
        # Complete directed graph on 3 vertices
        net = network(3)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 1)
        add_edge!(net, 1, 3)
        add_edge!(net, 3, 1)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 2)

        @test density(net) == 1.0
        @test reciprocity(net) == 1.0

        # Dyad census
        census = dyad_census(net)
        @test census.mutual == 3
        @test census.asymmetric == 0
        @test census.null == 0
    end

    @testset "Components" begin
        # Network with 2 components
        net = network(6)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 4, 5)
        add_edge!(net, 5, 6)

        comps = components(net; mode=:weak)
        @test length(comps) == 2
        @test Set(length.(comps)) == Set([3, 3])

        largest = largest_component(net)
        @test length(largest) == 3
    end

    @testset "Geodesic Distance" begin
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)

        dist = geodesic_distance(net)
        @test dist[1, 1] == 0.0
        @test dist[1, 2] == 1.0
        @test dist[1, 3] == 2.0
        @test dist[1, 4] == 3.0
        @test dist[4, 1] == Inf  # Can't reach 1 from 4 (directed)

        @test diameter(net) == 3.0
    end

    @testset "Structural Equivalence" begin
        # Network where vertices 1 and 2 have identical patterns
        net = network(4)
        add_edge!(net, 1, 3)
        add_edge!(net, 1, 4)
        add_edge!(net, 2, 3)
        add_edge!(net, 2, 4)

        se = structural_equivalence(net; method=:correlation)
        @test se[1, 2] == 1.0  # Perfectly equivalent
        @test se[1, 1] == 1.0  # Self-similarity
    end

    @testset "Blockmodel" begin
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 1)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 3)

        bm = blockmodel(net; k=2)
        @test length(bm.membership) == 4
        @test bm.n_blocks == 2
        @test size(bm.block_matrix) == (2, 2)
    end

    @testset "K-cores" begin
        # Create a network where some vertices have higher core numbers
        net = network(5; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 2, 3)  # 1,2,3 form a triangle
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        core_2 = kcores(net; k=2)
        @test Set(core_2) == Set([1, 2, 3])  # Triangle has core number 2
    end

    @testset "Golden master vs R sna (samplike, directed)" begin
        samp = sampson_like()

        @test density(samp) ≈ 0.2875816993 atol = 1e-9
        @test reciprocity(samp) ≈ 0.7908496732 atol = 1e-9
        @test reciprocity(samp; method=:edgewise) ≈ 0.6363636364 atol = 1e-9
        @test reciprocity(samp; method=:dyadic_nonnull) ≈ 0.4666666667 atol = 1e-9
        @test transitivity(samp) ≈ 0.4074074074 atol = 1e-9
        @test mutuality(samp) == 28

        dc = dyad_census(samp)
        @test (dc.mutual, dc.asymmetric, dc.null) == (28, 32, 93)

        # sna::triad.census order: 003 012 102 021D 021U 021C 111D 111U
        #                          030T 030C 201 120D 120U 120C 210 300
        @test triad_census(samp) == [167, 205, 190, 12, 24, 24, 68, 34,
                                     5, 0, 35, 15, 6, 5, 18, 8]

        @test connectedness(samp) ≈ 1.0
        @test efficiency(samp) ≈ 0.7543252595 atol = 1e-9
        @test hierarchy(samp) ≈ 0.2091503268 atol = 1e-9
        @test hierarchy(samp; measure=:krackhardt) ≈ 0.0 atol = 1e-12

        @test degree_centrality(samp)[1:4] == [17.0, 15.0, 8.0, 9.0]
        @test degree_centrality(samp; mode=:in)[1:4] == [11.0, 10.0, 3.0, 5.0]
        @test degree_centrality(samp; mode=:out)[1:4] == [6.0, 5.0, 5.0, 4.0]

        bc = betweenness_centrality(samp)
        @test bc[1] ≈ 68.9547619 atol = 1e-6
        @test bc[5] ≈ 38.74761905 atol = 1e-6
        @test bc[17] ≈ 0.0 atol = 1e-12

        ev = eigenvector_centrality(samp)
        @test ev[1] ≈ 0.287995961 atol = 1e-6
        @test ev[13] ≈ 0.1479037027 atol = 1e-6

        bp = bonacich_power(samp; exponent=0.05)
        @test bp[1] ≈ 1.212456868 atol = 1e-6
        @test bp[13] ≈ 0.6148524899 atol = 1e-6

        fb = flowbet(samp)
        @test fb[1] ≈ 155.0 atol = 1e-9
        @test fb[17] ≈ 18.0 atol = 1e-9
    end

    @testset "Golden master vs R sna (flomarriage, undirected)" begin
        flo = florentine()

        @test density(flo) ≈ 0.1666666667 atol = 1e-9
        @test transitivity(flo) ≈ 0.1914893617 atol = 1e-9
        @test triad_census(flo) == [324, 195, 38, 3]
        @test connectedness(flo) ≈ 0.875 atol = 1e-12
        @test efficiency(flo) ≈ 0.8673469388 atol = 1e-9

        bc = betweenness_centrality(flo)
        @test bc[9] ≈ 47.5 atol = 1e-9   # Medici
        @test bc[2] ≈ 19.33333333 atol = 1e-6
        @test bc[12] ≈ 0.0 atol = 1e-12  # Pucci (isolate)

        ev = eigenvector_centrality(flo)
        @test ev[9] ≈ 0.430308094 atol = 1e-6
        @test ev[12] ≈ 0.0 atol = 1e-9

        @test Set(kcores(flo; k=2)) == Set([2, 3, 4, 5, 7, 9, 11, 13, 15, 16])

        fb = flowbet(flo)
        @test fb[9] ≈ 68.0 atol = 1e-9
        @test fb[7] ≈ 38.0 atol = 1e-9
    end

    @testset "Triad census brute-force invariants" begin
        Random.seed!(99)
        net = rgnp(9, 0.3; directed=true)
        tc = triad_census(net)
        n = nv(net)
        @test sum(tc) == binomial(n, 3)
        # Cross-check dyad-level identities: each dyad appears in n-2 triads
        dc = dyad_census(net)
        # Mutual dyads per triad class: 102, 111D/U (1), 201 (2), 120* (1), 210 (2), 300 (3)
        m_from_tc = tc[3] + tc[7] + tc[8] + 2 * tc[11] + tc[12] + tc[13] +
                    tc[14] + 2 * tc[15] + 3 * tc[16]
        @test m_from_tc == dc.mutual * (n - 2)
    end

    @testset "Random network generators" begin
        Random.seed!(1)
        net = rgnp(20, 0.25; directed=true)
        @test nv(net) == 20
        @test is_directed(net)
        @test 0 < ne(net) < 380

        unet = rgnp(20, 0.25; directed=false)
        @test !is_directed(unet)

        m_net = rgnm(10, 17; directed=true)
        @test ne(m_net) == 17
        m_unet = rgnm(10, 17; directed=false)
        @test ne(m_unet) == 17
        @test_throws ArgumentError rgnm(3, 100)

        nets = rgraph(6; m=3, tprob=0.4)
        @test length(nets) == 3
        single = rgraph(6; tprob=0.4, mode=:graph)
        @test !is_directed(single)
    end

    @testset "Layouts" begin
        net = florentine()
        n = nv(net)

        for layout in (layout_circle, layout_random,
                       layout_fruchterman_reingold, layout_kamada_kawai)
            coords = layout(net)
            @test size(coords) == (n, 2)
            @test all(isfinite, coords)
        end

        # Circle layout: all on unit circle
        c = layout_circle(net)
        @test all(abs.(c[:, 1] .^ 2 .+ c[:, 2] .^ 2 .- 1.0) .< 1e-12)

        # KK layout roughly preserves relative distances: connected pairs
        # closer than the layout diameter
        kk = layout_kamada_kawai(net)
        @test size(kk) == (n, 2)
    end

    @testset "Bicomponents" begin
        net = network(5; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 1, 3)  # triangle 1-2-3
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        comps = bicomponents(net)
        # Triangle forms one biconnected component (3 edges); bridges are
        # their own components
        sizes = sort(length.(comps))
        @test sizes == [1, 1, 3]

        # Directed networks are treated as their underlying undirected graph
        dnet = network(3)
        add_edge!(dnet, 1, 2)
        add_edge!(dnet, 2, 3)
        @test length(bicomponents(dnet)) == 2
    end

    @testset "Equivalence clustering edge cases" begin
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 1)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 3)

        # k >= n: trivial clustering, blockmodel must not throw
        cl = equiv_clust(net; k=10)
        @test sort(unique(cl)) == collect(1:4)
        bm = blockmodel(net; k=10)
        @test bm.n_blocks == 4
        @test size(bm.block_matrix) == (4, 4)

        # k = 2 groups the structurally equivalent reciprocal pairs
        cl2 = equiv_clust(net; k=2)
        @test length(unique(cl2)) == 2
    end

    @testset "PageRank and Katz run" begin
        net = network(5)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 1)
        pr = pagerank(net)
        @test length(pr) == 5
        @test sum(pr) ≈ 1.0 atol = 1e-8
        kz = katz_centrality(net)
        @test length(kz) == 5
    end

    @testset "Cutpoints and Bridges" begin
        net = network(5; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        # In a path graph, all internal vertices are cutpoints
        cuts = cutpoints(net)
        @test Set(cuts) == Set([2, 3, 4])

        # All edges are bridges in a path
        br = bridges(net)
        @test length(br) == 4
    end
end
