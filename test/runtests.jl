using SNA
using Network
using Test

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
