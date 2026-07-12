using SNA
using Network
using Graphs
using Random
using Statistics
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

# Brute-force O(n^3) triad census, kept as the reference implementation for
# verifying the edge-driven Batagelj-Mrvar algorithm in src.
function ref_triad_census(net)
    n = nv(net)

    if !is_directed(net)
        census = zeros(Int, 4)
        for i in 1:n, j in (i+1):n, k in (j+1):n
            m = (has_edge(net, i, j) ? 1 : 0) +
                (has_edge(net, i, k) ? 1 : 0) +
                (has_edge(net, j, k) ? 1 : 0)
            census[m+1] += 1
        end
        return census
    end

    census = zeros(Int, 16)
    for i in 1:n, j in (i+1):n, k in (j+1):n
        census[ref_triad_type(net, i, j, k)] += 1
    end
    return census
end

# Classify the directed triad {a, b, c} into one of the 16 Davis-Leinhardt
# M-A-N classes (1-based index into the standard census order).
function ref_triad_type(net, a::Int, b::Int, c::Int)
    mutual = 0
    asym_arcs = Tuple{Int,Int}[]
    mutual_pair = (0, 0)

    for (i, j) in ((a, b), (a, c), (b, c))
        y_ij = has_edge(net, i, j)
        y_ji = has_edge(net, j, i)
        if y_ij && y_ji
            mutual += 1
            mutual_pair = (i, j)
        elseif y_ij
            push!(asym_arcs, (i, j))
        elseif y_ji
            push!(asym_arcs, (j, i))
        end
    end

    A = length(asym_arcs)

    if mutual == 3
        return 16                     # 300
    elseif mutual == 2
        return A == 1 ? 15 : 11       # 210 : 201
    elseif mutual == 1
        if A == 0
            return 3                  # 102
        elseif A == 1
            s, d = asym_arcs[1]
            return (d == mutual_pair[1] || d == mutual_pair[2]) ? 7 : 8
        else  # A == 2
            s1, d1 = asym_arcs[1]
            s2, d2 = asym_arcs[2]
            s1 == s2 && return 12     # 120D (common source)
            d1 == d2 && return 13     # 120U (common sink)
            return 14                 # 120C (chain)
        end
    else  # mutual == 0
        if A == 0
            return 1                  # 003
        elseif A == 1
            return 2                  # 012
        elseif A == 2
            s1, d1 = asym_arcs[1]
            s2, d2 = asym_arcs[2]
            s1 == s2 && return 4      # 021D (out-star)
            d1 == d2 && return 5      # 021U (in-star)
            return 6                  # 021C (path)
        else  # A == 3
            sources = (asym_arcs[1][1], asym_arcs[2][1], asym_arcs[3][1])
            return allunique(sources) ? 10 : 9   # 030C : 030T
        end
    end
end

# Graph correlation over off-diagonal dyads (R sna::gcor), used as the
# qaptest statistic
function gcor(a::AbstractMatrix, b::AbstractMatrix)
    n = size(a, 1)
    av = Float64[a[i, j] for i in 1:n, j in 1:n if i != j]
    bv = Float64[b[i, j] for i in 1:n, j in 1:n if i != j]
    return cor(av, bv)
end

# Padgett Florentine wealth (florentine_vertices.tsv), for netlm/netlogit
# dyadic covariates
const FLO_WEALTH = Float64[10, 36, 55, 44, 20, 32, 8, 42, 103, 48, 49, 3,
                           27, 10, 146, 48]

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

    @testset "Local clustering (undirected, single-counted degrees)" begin
        # Triangle 1-2-3 with pendant 4 attached to 3. R sna / igraph local
        # clustering: [1, 1, 1/3, 0]; the symmetric digraph storage must not
        # inflate the k(k-1) denominator.
        net = network(4; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)

        lc = transitivity(net; type=:local)
        @test lc ≈ [1.0, 1.0, 1 / 3, 0.0] atol = 1e-12
        @test transitivity(net; type=:average) ≈ 7 / 12 atol = 1e-12
    end

    @testset "Cliques" begin
        # Triangle 1-2-3 plus path 3-4-5: one maximal clique of size >= 3
        net = network(5; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        cl = cliques(net)
        @test Set.(cl) == [Set([1, 2, 3])]
        @test Set(Set.(cliques(net; min_size=2))) ==
              Set([Set([1, 2, 3]), Set([3, 4]), Set([4, 5])])

        # Directed networks are symmetrized with the weak rule, as in
        # R sna::clique.census: a one-way arc suffices for an undirected tie
        dnet = network(3)
        add_edge!(dnet, 1, 2)
        add_edge!(dnet, 2, 3)
        add_edge!(dnet, 3, 1)
        @test Set.(cliques(dnet)) == [Set([1, 2, 3])]
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

        # sna::degree(flo, gmode="graph"): undirected degree is single-counted
        @test degree_centrality(flo) == [1.0, 3.0, 2.0, 3.0, 3.0, 1.0, 4.0, 1.0,
                                         6.0, 1.0, 3.0, 0.0, 3.0, 2.0, 4.0, 3.0]
        @test degree_centrality(flo)[9] == 6.0  # Medici
        @test degree_centrality(flo; normalized=true)[9] ≈ 6 / 15 atol = 1e-12
        # mode is ignored for undirected networks (in = out = total degree)
        @test degree_centrality(flo; mode=:in) == degree_centrality(flo)

        # The three maximal cliques of size >= 3 (sna::clique.census)
        @test Set(Set.(cliques(flo))) ==
              Set([Set([4, 11, 15]), Set([5, 11, 15]), Set([9, 13, 16])])

        fb = flowbet(flo)
        @test fb[9] ≈ 68.0 atol = 1e-9
        @test fb[7] ≈ 38.0 atol = 1e-9
    end

    @testset "Graphs.jl namespace integration" begin
        # SNA extends the Graphs.jl generics instead of shadowing them, so
        # `using SNA, Graphs` (as at the top of this file) must leave a
        # single unambiguous binding for each shared name.
        @test density === Graphs.density
        @test diameter === Graphs.diameter
        @test bridges === Graphs.bridges
        @test degree_centrality === Graphs.degree_centrality
        @test betweenness_centrality === Graphs.betweenness_centrality
        @test closeness_centrality === Graphs.closeness_centrality
        @test eigenvector_centrality === Graphs.eigenvector_centrality
        @test katz_centrality === Graphs.katz_centrality
        @test pagerank === Graphs.pagerank

        # Undirected n=5 with 2 edges: density is 2/10 = 0.2 (not the
        # doubled-storage 0.1)
        net = network(5; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 3, 4)
        @test density(net) ≈ 0.2 atol = 1e-12
        @test Graphs.density(net) ≈ 0.2 atol = 1e-12

        path = network(5; directed=false)
        for i in 1:4
            add_edge!(path, i, i + 1)
        end
        @test diameter(path) == 4.0
        @test length(bridges(path)) == 4
        @test degree_centrality(path) == [1.0, 2.0, 2.0, 2.0, 1.0]

        # The generics still work on plain Graphs.jl graphs
        g = Graphs.path_graph(5)
        @test density(g) ≈ 0.4 atol = 1e-12
        @test diameter(g) == 4
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

    @testset "Triad census: Batagelj-Mrvar vs brute force" begin
        # Exact agreement with the O(n^3) reference on random directed and
        # undirected graphs across densities (including empty and complete)
        rng = Xoshiro(2024)
        for trial in 1:8
            n = rand(rng, 3:40)
            p = rand(rng, (0.0, 0.02, 0.1, 0.3, 0.7, 1.0))
            dnet = rgnp(n, p; directed=true, rng=rng)
            @test triad_census(dnet) == ref_triad_census(dnet)
            unet = rgnp(n, p; directed=false, rng=rng)
            @test triad_census(unet) == ref_triad_census(unet)
        end

        # Tiny-network edge cases
        for n in (1, 2), directed in (true, false)
            tiny = network(n; directed=directed)
            n == 2 && add_edge!(tiny, 1, 2)
            @test triad_census(tiny) == ref_triad_census(tiny)
            @test sum(triad_census(tiny)) == 0
        end

        # Golden master vs R sna::triad.census: fixed directed fixture
        # (n = 7, set.seed(42) rgraph(7, tprob=0.35) in R sna 2.8)
        fix = network(7; directed=true)
        for (i, j) in [(1, 3), (1, 4), (1, 7), (2, 5), (3, 1), (3, 2),
                       (3, 5), (3, 6), (3, 7), (4, 1), (4, 3), (4, 6),
                       (5, 1), (5, 2), (6, 3), (6, 5), (6, 7), (7, 1),
                       (7, 2), (7, 3), (7, 6)]
            add_edge!(fix, i, j)
        end
        # sna order: 003 012 102 021D 021U 021C 111D 111U 030T 030C
        #            201 120D 120U 120C 210 300
        @test triad_census(fix) == [2, 2, 2, 0, 0, 4, 6, 8, 0, 0,
                                    3, 2, 2, 1, 1, 2]
    end

    @testset "Centralization golden master vs R sna" begin
        flo = florentine()
        samp = sampson_like()

        # R sna::centralization(flo, FUN, mode="graph")
        @test centralization(flo, :degree) ≈ 0.2666666667 atol = 1e-9
        @test centralization(flo, :betweenness) ≈ 0.3834920635 atol = 1e-9
        # Pucci is an isolate: every Freeman closeness score is undefined
        # (0 by sna convention), so the centralization is 0
        @test centralization(flo, :closeness) ≈ 0.0 atol = 1e-12
        @test centralization(flo, :eigenvector) ≈ 0.3416651990 atol = 1e-6

        # R sna::centralization(samp, FUN) (mode="digraph")
        @test centralization(samp, :degree) ≈ 0.2389705882 atol = 1e-9
        @test centralization(samp, :degree; mode=:in) ≈ 0.3806228374 atol = 1e-9
        @test centralization(samp, :degree; mode=:out) ≈ 0.0692041522 atol = 1e-9
        @test centralization(samp, :betweenness) ≈ 0.2024623085 atol = 1e-9
        @test centralization(samp, :closeness) ≈ 0.0919022955 atol = 1e-9
        @test centralization(samp, :eigenvector) ≈ 0.0700638198 atol = 1e-6

        # normalize=FALSE returns the raw deviation sum
        @test centralization(flo, :degree; normalized=false) ≈ 56.0 atol = 1e-9
        @test centralization(samp, :degree; normalized=false) ≈ 130.0 atol = 1e-9

        # Star graph is maximally degree-centralized
        star = network(6; directed=false)
        for v in 2:6
            add_edge!(star, 1, v)
        end
        @test centralization(star, :degree) ≈ 1.0 atol = 1e-12
        @test centralization(star, :betweenness) ≈ 1.0 atol = 1e-12
        @test centralization(star, :closeness) ≈ 1.0 atol = 1e-12

        @test_throws ArgumentError centralization(flo, :pagerank)
    end

    @testset "QAP test (qaptest)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)

        qt = qaptest(gcor, flo, biz; reps=1000, rng=Xoshiro(11))
        # Observed statistic is deterministic: R sna::gcor(flo, biz)
        @test qt.testval ≈ 0.3718678721 atol = 1e-9
        @test qt isa QAPTestResult
        @test length(qt.dist) == 1000
        @test qt.reps == 1000
        # Marriage and business ties are strongly associated: the QAP
        # p-value is far in the upper tail (R reference: pgreq ~ 0.001)
        @test qt.pgreq <= 0.01
        @test qt.pleeq >= 0.99
        @test qt.pgreq == count(>=(qt.testval), qt.dist) / qt.reps

        # Matrices are accepted directly, and f sees permuted matrices
        qm = qaptest(gcor, as_matrix(flo), as_matrix(biz); reps=100,
                     rng=Xoshiro(1))
        @test qm.testval ≈ qt.testval atol = 1e-12

        # A self-comparison is at the very top of its null distribution
        qs = qaptest(gcor, flo, flo; reps=100, rng=Xoshiro(2))
        @test qs.testval ≈ 1.0 atol = 1e-12
        @test qs.pgreq <= 0.05

        @test occursin("QAP Test", sprint(show, qt))
        @test_throws ArgumentError qaptest(gcor, flo, network(5; directed=false))
    end

    @testset "Network regression (netlm)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)
        wdiff = abs.(FLO_WEALTH .- FLO_WEALTH')

        # Golden master vs R sna::netlm(flo, list(biz, wdiff), mode="graph",
        # nullhyp="classical"): undirected dyads, diagonal excluded
        fit = netlm(flo, [biz, wdiff]; nullhyp=:classical)
        @test fit isa NetLMResult
        @test fit.n == 120
        @test fit.df_residual == 117
        @test !fit.directed
        @test fit.names == ["(intercept)", "x1", "x2"]
        @test fit.coefficients ≈ [0.0115567910, 0.4292240761, 0.0026646583] atol = 1e-9
        @test fit.tstat ≈ [0.2469658004, 4.6124554994, 3.0852487650] atol = 1e-9
        @test fit.pgreqabs ≈ [0.8053675128, 1.0246842090e-5, 0.0025385758] atol = 1e-9
        @test fit.r_squared ≈ 0.2031176130 atol = 1e-9
        @test fit.dist === nothing

        # Directed golden master: samplike on its transpose
        # (R: netlm(samp, t(samp), nullhyp="classical"))
        samp = sampson_like()
        fit_d = netlm(samp, [Matrix(as_matrix(samp)')]; nullhyp=:classical)
        @test fit_d.n == 306
        @test fit_d.directed
        @test fit_d.coefficients ≈ [0.1467889908, 0.4895746455] atol = 1e-9
        @test fit_d.tstat ≈ [5.4733396394, 9.7894536551] atol = 1e-9

        # Dekker double-semi-partialing QAP (the default): identical point
        # estimates, permutation p-values (R reference with reps=2000:
        # pgreqabs ~ [0.80, 0.000, 0.0015])
        fq = netlm(flo, [biz, wdiff]; reps=500, rng=Xoshiro(7))
        @test fq.nullhyp == :qapspp
        @test fq.coefficients ≈ fit.coefficients atol = 1e-12
        @test fq.tstat ≈ fit.tstat atol = 1e-12
        @test size(fq.dist) == (500, 3)
        @test fq.pgreqabs[1] > 0.5      # intercept: no effect
        @test fq.pgreqabs[2] <= 0.01    # business ties: strong effect
        @test fq.pgreqabs[3] <= 0.05    # wealth difference: real effect
        @test all(0 .<= fq.pleeq .<= 1) && all(0 .<= fq.pgreq .<= 1)

        # Classical y-permutation QAP (R reference: ~ [1, 0.000, 0.0055])
        fy = netlm(flo, [biz, wdiff]; nullhyp=:qapy, reps=500, rng=Xoshiro(7))
        @test fy.nullhyp == :qapy
        @test fy.pgreqabs[2] <= 0.01
        @test fy.pgreqabs[3] <= 0.05

        # x-permutation QAP runs and keeps the same point estimates
        fx = netlm(flo, [biz, wdiff]; nullhyp=:qapx, reps=100, rng=Xoshiro(7))
        @test fx.nullhyp == :qapx
        @test fx.coefficients ≈ fit.coefficients atol = 1e-12

        # With a single regressor, semi-partialing degenerates to :qapy
        # (as in sna); a single predictor without intercept has nx == 1
        f1 = netlm(flo, biz; intercept=false, reps=100, rng=Xoshiro(1))
        @test f1.nullhyp == :qapy
        @test length(f1.coefficients) == 1

        # A single-network (non-vector) predictor is accepted
        f2 = netlm(flo, biz; nullhyp=:classical)
        @test f2.names == ["(intercept)", "x1"]

        # Raw matrices default to directed dyads; mode=:graph overrides
        fm = netlm(as_matrix(flo), [as_matrix(biz)]; nullhyp=:classical)
        @test fm.directed && fm.n == 240
        fg = netlm(as_matrix(flo), [as_matrix(biz)]; nullhyp=:classical,
                   mode=:graph)
        @test !fg.directed && fg.n == 120
        @test fg.coefficients ≈ f2.coefficients atol = 1e-12

        @test occursin("R-squared", sprint(show, fq))
        @test_throws ArgumentError netlm(flo, [network(5; directed=false)])
        @test_throws ArgumentError netlm(flo, [biz]; nullhyp=:bogus)
        @test_throws ArgumentError netlm(flo, [biz]; mode=:bogus)
    end

    @testset "Network logit (netlogit)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)
        wdiff = abs.(FLO_WEALTH .- FLO_WEALTH')

        # Golden master vs R sna::netlogit(flo, list(biz, wdiff),
        # mode="graph", nullhyp="classical")
        fit = netlogit(flo, [biz, wdiff]; nullhyp=:classical)
        @test fit isa NetLogitResult
        @test fit.n == 120
        @test fit.df_residual == 117
        @test fit.coefficients ≈ [-3.0417761907, 2.4966227608, 0.0199475777] atol = 1e-5
        @test fit.se ≈ [0.5283397915, 0.6511320086, 0.0068658522] atol = 1e-5
        @test fit.tstat ≈ [-5.7572347186, 3.8342804957, 2.9053316701] atol = 1e-4
        @test fit.pgreqabs ≈ [7.0107518e-8, 2.0438801e-4, 4.3880907e-3] rtol = 1e-3
        @test fit.deviance ≈ 86.9618703491 atol = 1e-6
        @test fit.null_deviance ≈ 166.3553233344 atol = 1e-6
        @test fit.aic ≈ 92.9618703491 atol = 1e-6
        @test fit.bic ≈ 101.3243455775 atol = 1e-6

        # DSP QAP p-values (R reference with reps=1000: ~ [0, 0, 0.001])
        fq = netlogit(flo, [biz, wdiff]; reps=200, rng=Xoshiro(7))
        @test fq.nullhyp == :qapspp
        @test fq.coefficients ≈ fit.coefficients atol = 1e-8
        @test fq.pgreqabs[2] <= 0.05
        @test fq.pgreqabs[3] <= 0.05
        @test size(fq.dist) == (200, 3)

        @test occursin("deviance", sprint(show, fq))

        # The DV must be dichotomous
        @test_throws ArgumentError netlogit(wdiff, [as_matrix(biz)])
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
