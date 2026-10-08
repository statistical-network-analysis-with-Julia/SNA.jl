using SNA
using NetworkCore
using Graphs
using Random
using Statistics
using LinearAlgebra
using Test
using Aqua

const R_SNA = NetworkCore.load_golden(joinpath(@__DIR__, "fixtures", "sna_reference.toml"))
const R_FUZZ = NetworkCore.load_golden(joinpath(@__DIR__, "fixtures", "sna_fuzz.toml"))
const R_INF = NetworkCore.load_golden(joinpath(@__DIR__, "fixtures", "sna_inference.toml"))
sampson_like() = network_from_matrix(reshape(R_SNA.values["samp_adjacency"], 18, 18); directed=true)
florentine() = network_from_matrix(reshape(R_SNA.values["flo_adjacency"], 16, 16); directed=false)

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

# Padgett Florentine wealth (florentine_vertices.tsv), for netlm/netlogit
# dyadic covariates
const FLO_WEALTH = Float64.(R_SNA.values["flo_wealth"])

@testset "SNA.jl" begin
    @testset "Degree Centrality" begin
        net = network(5)
        add_edge!(net, 1, 2)
        add_edge!(net, 1, 3)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        dc = degreecent(net; cmode=:outdegree)
        @test dc[1] == 2.0  # 1 sends to 2, 3
        @test dc[3] == 1.0  # 3 sends to 4
        @test dc[5] == 0.0  # 5 sends to nobody

        dc_in = degreecent(net; cmode=:indegree)
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

        @test gden(net) == 1.0
        @test grecip(net) == 1.0

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

        cd = component_dist(net; connected=:weak)
        @test cd.membership == [1, 1, 1, 2, 2, 2]
        @test cd.csize == [3, 3]
        @test cd.cdist == [0, 0, 2, 0, 0, 0]
        # The development-time name `components` (it collided with
        # Distributions.components) was never released and has no alias
        @test !isdefined(SNA, :components)

        largest = largest_component(net; connected=:weak)
        @test length(largest) == 3
        @test length(largest_component(net)) == 1
    end

    @testset "Geodesic Distance" begin
        net = network(4)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)

        dist = geodist(net).gdist
        @test dist[1, 1] == 0.0
        @test dist[1, 2] == 1.0
        @test dist[1, 3] == 2.0
        @test dist[1, 4] == 3.0
        @test dist[4, 1] == Inf  # Can't reach 1 from 4 (directed)

        @test maximum(dist) == Inf # directed endpoints cannot reach each other
        @test geodist(net).counts[1, 4] == 1
        @test geodist(net; inf_replace=0).gdist[4, 1] == 0
        @test !isdefined(SNA, :geodesic_distance)    # renamed geodist, no alias
    end

    @testset "Structural Equivalence" begin
        # Network where vertices 1 and 2 have identical patterns
        net = network(4)
        add_edge!(net, 1, 3)
        add_edge!(net, 1, 4)
        add_edge!(net, 2, 3)
        add_edge!(net, 2, 4)

        se = sedist(net; method=:correlation)
        @test se[1, 2] == 1.0  # Perfectly equivalent
        @test se[1, 1] == 1.0  # Self-similarity
        @test sedist(net)[1, 2] == 0.0          # Hamming (sna's default) is a count

        # The profile of a pair leaves out the self cells and
        # the ties between the two actors (positions i, j, n+i, n+j), as
        # sna::sedist does. Undirected triangle 1-2-3 plus isolate 4: triangle
        # members are equivalent; directed: 1 and 2 are tied to each other and
        # identical towards 3 and 4. Values are R sna 2.8's.
        tri = network(4; directed=false)
        add_edges!(tri, [(1, 2), (1, 3), (2, 3)])
        pair = network(4)
        add_edges!(pair, [(1, 2), (2, 1), (1, 3), (2, 3), (4, 1), (4, 2)])
        for (id, g) in (("triangle", tri), ("pair", pair)),
            m in (:hamming, :correlation, :euclidean)
            @test NetworkCore.check_golden(R_FUZZ, "se_$(id)_$m", vec(sedist(g; method=m)))
        end
        @test sedist(tri; method=:euclidean)[1, 2] == 0.0
        @test sedist(tri; method=:euclidean)[1, 4] == 2.0
        @test sedist(pair; method=:correlation)[1, 2] == 1.0
        # diag=true keeps every cell, the pre-0.2 profile
        @test sedist(pair; diag=true)[1, 2] == 4.0
        # ...and the clustering built on it now groups the equivalent pair
        @test equiv_clust(pair; k=3)[1] == equiv_clust(pair; k=3)[2]

        @test !isdefined(SNA, :structural_equivalence)   # renamed sedist, no alias
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

        @test kcores(net) == [2, 2, 2, 1, 1]          # core numbers, as R sna::kcores
        @test findall(>=(2), kcores(net)) == [1, 2, 3]
        # The pre-0.2 `k` keyword (k-core members) is refused with the new idiom
        err = try kcores(net; k=2) catch e e end
        @test err isa ArgumentError && occursin("findall(>=(2), kcores(net))", sprint(showerror, err))
        # Self-loops are ignored (sna diag=FALSE); Graphs.core_number would throw
        looped = network(5; directed=false, loops=true)
        add_edges!(looped, [(1, 2), (1, 3), (2, 3), (3, 4), (4, 5), (1, 1), (4, 4)])
        @test kcores(looped) == [2, 2, 2, 1, 1]
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

        # A vertex of degree < 2 has no coefficient (NaN), and the average
        # leaves it out, as igraph does
        lc = transitivity(net; type=:local)
        @test lc[1:3] ≈ [1.0, 1.0, 1 / 3] atol = 1e-12
        @test isnan(lc[4])
        @test transitivity(net; type=:average) ≈ 7 / 9 atol = 1e-12
        @test isnan(transitivity(network(3; directed=false); type=:average))
        # Self-loops are not neighbours
        looped = network(3; directed=false, loops=true)
        add_edges!(looped, [(1, 2), (1, 3), (2, 3), (1, 1)])
        @test transitivity(looped; type=:local) == [1.0, 1.0, 1.0]
    end

    @testset "Directed local clustering: Fagiolo (2007) total coefficient" begin
        # The complete digraph on three actors is fully
        # clustered (Graphs.jl's mixed definition gave 1/6)
        k3 = network(3)
        add_edges!(k3, [(1, 2), (2, 1), (1, 3), (3, 1), (2, 3), (3, 2)])
        @test transitivity(k3; type=:local) == [1.0, 1.0, 1.0]
        # Brute force: t_i = ½ Σ_{j,h} (a_ij+a_ji)(a_ih+a_hi)(a_jh+a_hj),
        # T_i = d_tot(d_tot − 1) − 2 d↔
        function fagiolo(A)
            n = size(A, 1)
            c = fill(NaN, n)
            for i in 1:n
                t = 0.0
                for j in 1:n, h in 1:n
                    (j == i || h == i || j == h) && continue
                    t += (A[i, j] + A[j, i]) * (A[i, h] + A[h, i]) * (A[j, h] + A[h, j])
                end
                dtot = sum(A[i, :]) + sum(A[:, i])
                dbi = sum(A[i, j] * A[j, i] for j in 1:n)
                T = dtot * (dtot - 1) - 2dbi
                T > 0 && (c[i] = t / 2 / T)
            end
            return c
        end
        rng = Xoshiro(2007)
        for _ in 1:40
            n = rand(rng, 3:9)
            g = rgnp(n, rand(rng, (0.2, 0.4, 0.7)); directed=true, rng=rng)
            add_edge!(g, 1, 2)                 # never empty
            A = Float64.(SNA._sociomatrix(g))
            @test isequal(round.(transitivity(g; type=:local); digits=12),
                          round.(fagiolo(A); digits=12))
        end
        # On a symmetric digraph it is the undirected Watts–Strogatz value
        sym = network(4)
        und = network(4; directed=false)
        for (i, j) in [(1, 2), (1, 3), (2, 3), (3, 4)]
            add_edges!(sym, [(i, j), (j, i)])
            add_edge!(und, i, j)
        end
        @test isequal(transitivity(sym; type=:local), transitivity(und; type=:local))
        # cmode=:weak is igraph's convention (direction ignored)
        @test transitivity(k3; type=:local, cmode=:weak) == [1.0, 1.0, 1.0]
        @test_throws ArgumentError transitivity(k3; type=:local, cmode=:bogus)
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
        @test isempty(cliques(dnet))
        @test Set.(cliques(dnet; symmetrize=:weak)) == [Set([1, 2, 3])]
    end

    @testset "Provenanced R sna reference: full vectors and directed cohesion" begin
        masks(sets) = sort([sum(2^(v-1) for v in unique(xs)) for xs in sets])
        for (id, net) in (("flo", florentine()), ("samp", sampson_like()))
            check(key, value) = NetworkCore.check_golden(R_SNA, id * "_" * key, value)
            @test check("density", gden(net))
            @test check("reciprocity", grecip(net))
            @test check("reciprocity_edgewise", grecip(net; measure=:edgewise))
            @test check("reciprocity_nonnull", grecip(net; measure=:dyadic_nonnull))
            @test check("transitivity", transitivity(net))
            @test check("dyad_census", collect(values(dyad_census(net))))
            @test check("triad_census", triad_census(net))
            @test mutuality(net) == dyad_census(net).mutual
            @test check("connectedness", connectedness(net))
            @test check("efficiency", efficiency(net))
            @test check("hierarchy", hierarchy(net))
            @test check("hierarchy_krackhardt", hierarchy(net; measure=:krackhardt))
            @test check("degree", degreecent(net))
            @test check("degree_in", degreecent(net; cmode=:indegree))
            @test check("degree_out", degreecent(net; cmode=:outdegree))
            @test check("betweenness", betweenness(net))
            @test check("closeness", closeness(net))
            @test check("eigenvector", evcent(net))
            @test check("bonacich", bonpow(net; exponent=0.05))
            @test check("flowbet", flowbet(net))
            @test check("strong_sizes", sort(component_dist(net).csize; rev=true))
            @test check("weak_sizes", sort(component_dist(net; connected=:weak).csize; rev=true))
            @test check("cutpoints", sort(cutpoints(net)))
            @test check("cutpoints_weak", sort(cutpoints(net; connected=:weak)))
            @test check("cutpoints_recursive", sort(cutpoints(net; connected=:recursive)))
            @test check("bicomponents", masks([collect(Iterators.flatten(c)) for c in bicomponents(net)]))
            @test check("cliques", masks(cliques(net)))
            for k in (1, 2, 5, 7)
                @test check("kcore_$k", findall(>=(k), kcores(net)))
            end
            for (measure, ref) in ((:degree,"degree"), (:betweenness,"betweenness"),
                                   (:closeness,"closeness"), (:eigenvector,"evcent"))
                @test check("centralization_" * ref, centralization(net, measure))
            end
            @test check("centralization_in", centralization(net, :degree; mode=:in))
            @test check("centralization_out", centralization(net, :degree; mode=:out))
            @test check("centralization_raw", centralization(net, :degree; normalized=false))
        end
        flo = florentine()
        # sna's closeness of a disconnected graph is 0 everywhere; its
        # harmonic variant is positive for every vertex except the isolate
        # (Pucci, vertex 12)
        @test NetworkCore.check_golden(R_SNA, "flo_closeness", closeness(flo))
        hc = closeness(flo; cmode=:suminvundir)
        @test findall(iszero, hc) == [12]
        cycle = network(4)
        add_edges!(cycle, [(1,2), (2,3), (3,1)])
        @test NetworkCore.check_golden(R_SNA, "cycle_cutpoints", cutpoints(cycle))
        @test NetworkCore.check_golden(R_SNA, "cycle_cutpoints_recursive", cutpoints(cycle; connected=:recursive))
        @test NetworkCore.check_golden(R_SNA, "cycle_strong_sizes", sort(component_dist(cycle).csize; rev=true))
        @test NetworkCore.check_golden(R_SNA, "cycle_closeness", closeness(cycle))
        @test NetworkCore.check_golden(R_SNA, "cycle_cliques", masks(cliques(cycle)))
    end

    @testset "No type piracy: Graphs.jl functions keep Graphs.jl semantics" begin
        # SNA adds no method to a Graphs.jl function. Its R-semantics measures
        # are SNA-owned names, and Graphs' own functions answer on a Network
        # exactly as they answer on the equivalent SimpleGraph/SimpleDiGraph.
        for f in (Graphs.degree_centrality, Graphs.betweenness_centrality,
                  Graphs.closeness_centrality, Graphs.eigenvector_centrality,
                  Graphs.katz_centrality, Graphs.pagerank, Graphs.density,
                  Graphs.diameter, Graphs.bridges)
            @test !any(m -> m.module === SNA, methods(f))
        end
        for (sna_name, graphs_name) in ((:degreecent, :degree_centrality),
                                        (:betweenness, :betweenness_centrality),
                                        (:closeness, :closeness_centrality),
                                        (:evcent, :eigenvector_centrality),
                                        (:gden, :density))
            @test getfield(SNA, sna_name) !== getfield(Graphs, graphs_name)
            @test !Base.isexported(SNA, graphs_name)
        end
        # A disconnected undirected network: Graphs' normalized betweenness and
        # component-scaled closeness on a Network equal Graphs' own answers on
        # the SimpleGraph, while SNA's measures follow R sna.
        net = network(5; directed=false)
        add_edges!(net, [(1, 2), (2, 3), (4, 5)])
        sg = Graphs.SimpleGraph(5)
        for (i, j) in [(1, 2), (2, 3), (4, 5)]
            Graphs.add_edge!(sg, i, j)
        end
        @test Graphs.betweenness_centrality(net) ≈ Graphs.betweenness_centrality(sg)
        @test Graphs.closeness_centrality(net) ≈ Graphs.closeness_centrality(sg)
        @test Graphs.degree_centrality(net) ≈ Graphs.degree_centrality(sg)
        @test Graphs.betweenness_centrality(net; normalize=false) ≈ [0, 1, 0, 0, 0]
        @test betweenness(net) == [0, 1, 0, 0, 0]
        @test all(iszero, closeness(net))                 # sna: unreachable → 0
        @test gden(net) ≈ 0.3 atol = 1e-12
        @test Graphs.density(net) ≈ 0.3 atol = 1e-12

        # Co-loading Graphs, Distributions, NetworkCore and SNA leaves every name
        # each of them exports defined (no conflicting exports)
        code = """
            using Graphs, Distributions, NetworkCore, SNA
            bad = [(m, n) for m in (Graphs, Distributions, NetworkCore, SNA)
                   for n in names(m) if Base.isexported(m, n) && !isdefined(Main, n)]
            isempty(bad) || error("names left undefined by co-loading: \$bad")
            print("ok")
            """
        cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`
        @test read(cmd, String) == "ok"
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

    end

    @testset "Centralization golden master vs R sna" begin
        flo = florentine()
        samp = sampson_like()

        # Star graph is maximally degree-centralized
        star = network(6; directed=false)
        for v in 2:6
            add_edge!(star, 1, v)
        end
        @test centralization(star, :degree) ≈ 1.0 atol = 1e-12
        @test centralization(star, :betweenness) ≈ 1.0 atol = 1e-12
        @test centralization(star, :closeness) ≈ 1.0 atol = 1e-12

        @test_throws ArgumentError centralization(flo, :pagerank)
        # The function form passes keywords to the measure, as R passes `...`
        @test centralization(flo, degreecent) ≈ centralization(flo, :degree)
        @test centralization(samp, degreecent; cmode=:indegree) ≈
              centralization(samp, :degree; mode=:in)
        @test isnan(centralization(network(2; directed=false), degreecent))   # 0/0, as sna
        @test_throws ArgumentError centralization(flo, infocent)
        @test centralization(flo, infocent; normalized=false) ≥ 0
    end

    @testset "QAP test (qaptest)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)

        qt = qaptest(gcor, flo, biz; n_sim=1000, rng=Xoshiro(11))
        # Observed statistic is deterministic: R sna::gcor(flo, biz)
        @test NetworkCore.check_golden(R_SNA, "qap_gcor", qt.testval)
        @test qt isa QAPTestResult
        @test length(qt.dist) == 1000
        @test qt.reps == 1000
        # Marriage and business ties are strongly associated: the QAP
        # p-value is far in the upper tail (R reference: pgreq ~ 0.001)
        @test qt.pgreq <= 0.01
        @test qt.pleeq >= 0.99
        @test qt.pgreq == count(>=(qt.testval), qt.dist) / qt.reps

        # Matrices are accepted directly, and f sees permuted matrices
        qm = qaptest(gcor, as_matrix(flo), as_matrix(biz); n_sim=100,
                     rng=Xoshiro(1))
        @test qm.testval ≈ qt.testval atol = 1e-12

        # A self-comparison is at the very top of its null distribution
        qs = qaptest(gcor, flo, flo; n_sim=100, rng=Xoshiro(2))
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
        @test NetworkCore.check_golden(R_SNA, "lm_coefficients", fit.coefficients)
        @test NetworkCore.check_golden(R_SNA, "lm_tstat", fit.tstat)
        @test NetworkCore.check_golden(R_SNA, "lm_pgreqabs", fit.pgreqabs)
        @test NetworkCore.check_golden(R_SNA, "lm_r_squared", fit.r_squared)
        @test fit.dist === nothing

        # Directed golden master: samplike on its transpose
        # (R: netlm(samp, t(samp), nullhyp="classical"))
        samp = sampson_like()
        fit_d = netlm(samp, [Matrix(as_matrix(samp)')]; nullhyp=:classical)
        @test fit_d.n == 306
        @test fit_d.directed
        @test NetworkCore.check_golden(R_SNA, "lm_directed_coefficients", fit_d.coefficients)
        @test NetworkCore.check_golden(R_SNA, "lm_directed_tstat", fit_d.tstat)

        # Dekker double-semi-partialing QAP (the default): identical point
        # estimates, permutation p-values (R reference with reps=2000:
        # pgreqabs ~ [0.80, 0.000, 0.0015])
        fq = netlm(flo, [biz, wdiff]; n_sim=500, rng=Xoshiro(7))
        @test fq.nullhyp == :qapspp
        @test fq.coefficients ≈ fit.coefficients atol = 1e-12
        @test fq.tstat ≈ fit.tstat atol = 1e-12
        @test size(fq.dist) == (500, 3)
        @test fq.pgreqabs[1] > 0.5      # intercept: no effect
        @test fq.pgreqabs[2] <= 0.01    # business ties: strong effect
        @test fq.pgreqabs[3] <= 0.05    # wealth difference: real effect
        @test all(0 .<= fq.pleeq .<= 1) && all(0 .<= fq.pgreq .<= 1)

        # Classical y-permutation QAP (R reference: ~ [1, 0.000, 0.0055])
        fy = netlm(flo, [biz, wdiff]; nullhyp=:qapy, n_sim=500, rng=Xoshiro(7))
        @test fy.nullhyp == :qapy
        @test fy.pgreqabs[2] <= 0.01
        @test fy.pgreqabs[3] <= 0.05

        # x-permutation QAP runs and keeps the same point estimates
        fx = netlm(flo, [biz, wdiff]; nullhyp=:qapx, n_sim=100, rng=Xoshiro(7))
        @test fx.nullhyp == :qapx
        @test fx.coefficients ≈ fit.coefficients atol = 1e-12

        # With a single regressor, semi-partialing degenerates to :qapy
        # (as in sna); a single predictor without intercept has nx == 1
        f1 = netlm(flo, biz; intercept=false, n_sim=100, rng=Xoshiro(1))
        @test f1.nullhyp == :qapy
        @test length(f1.coefficients) == 1

        # A single-network (non-vector) predictor is accepted
        f2 = netlm(flo, biz; nullhyp=:classical)
        @test f2.names == ["(intercept)", "x1"]

        # Raw matrices: mode=:auto treats symmetric data as undirected (stats
        # S7: sna's "digraph" default counts every pair twice and makes the
        # classical standard errors √2 too small); mode=:digraph is sna's
        fm = netlm(as_matrix(flo), [as_matrix(biz)]; nullhyp=:classical)
        @test !fm.directed && fm.n == 120
        @test fm.coefficients ≈ f2.coefficients atol = 1e-12
        @test fm.tstat ≈ f2.tstat atol = 1e-12
        fd = netlm(as_matrix(flo), [as_matrix(biz)]; nullhyp=:classical,
                   mode=:digraph)
        @test fd.directed && fd.n == 240
        @test fd.coefficients ≈ f2.coefficients atol = 1e-12
        @test fd.tstat ≈ sqrt(238 / 118) .* f2.tstat atol = 1e-9     # the √2 of S7
        # An asymmetric predictor makes the data directed
        asym = Matrix{Float64}(as_matrix(biz)); asym[1, 2] = 1 - asym[1, 2]
        @test netlm(as_matrix(flo), asym; nullhyp=:classical).directed
        # R's keyword name `reps` is `n_sim` here; there is no alias
        @test_throws MethodError netlm(flo, biz; reps=20, rng=Xoshiro(1))

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
        @test NetworkCore.check_golden(R_SNA, "logit_coefficients", fit.coefficients)
        @test NetworkCore.check_golden(R_SNA, "logit_se", fit.se)
        @test NetworkCore.check_golden(R_SNA, "logit_tstat", fit.tstat)
        @test NetworkCore.check_golden(R_SNA, "logit_pgreqabs", fit.pgreqabs)
        @test NetworkCore.check_golden(R_SNA, "logit_deviance", fit.deviance)
        @test NetworkCore.check_golden(R_SNA, "logit_null_deviance", fit.null_deviance)
        @test NetworkCore.check_golden(R_SNA, "logit_aic", fit.aic)
        @test NetworkCore.check_golden(R_SNA, "logit_bic", fit.bic)

        # DSP QAP p-values (R reference with reps=1000: ~ [0, 0, 0.001])
        fq = netlogit(flo, [biz, wdiff]; n_sim=200, rng=Xoshiro(7))
        @test fq.nullhyp == :qapspp
        @test fq.statistic == :lr
        @test fq.coefficients ≈ fit.coefficients atol = 1e-8
        @test fq.pgreqabs[2] <= 0.05
        @test fq.pgreqabs[3] <= 0.05
        @test size(fq.dist) == (200, 3)
        @test occursin("LR z", sprint(show, fq))
        # The LR statistic is the signed root of the deviance drop
        idx = SNA._dyad_indices(16, false)
        X = hcat(ones(120), SNA._gvectorize(SNA._sociomatrix(biz), idx), SNA._gvectorize(wdiff, idx))
        yv = SNA._gvectorize(SNA._sociomatrix(flo), idx)
        for k in 1:3
            d_red = SNA._logit_fit(X[:, setdiff(1:3, k)], yv).deviance
            @test fq.tstat[k] ≈ sign(fit.coefficients[k]) * sqrt(d_red - fit.deviance) atol = 1e-6
        end
        # sna's Wald statistic is still available
        fw = netlogit(flo, [biz, wdiff]; statistic=:wald, n_sim=100, rng=Xoshiro(7))
        @test fw.tstat ≈ fit.tstat atol = 1e-6
        @test_throws ArgumentError netlogit(flo, biz; statistic=:bogus)

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
        @test sizes == [3]
        @test sort(length.(bicomponents(net; min_size=2))) == [1, 1, 3]

        # Directed networks are treated as their underlying undirected graph
        dnet = network(3)
        add_edge!(dnet, 1, 2)
        add_edge!(dnet, 2, 3)
        @test isempty(bicomponents(dnet))
        @test length(bicomponents(dnet; symmetrize=:weak, min_size=2)) == 2
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

    @testset "Cutpoints" begin
        net = network(5; directed=false)
        add_edge!(net, 1, 2)
        add_edge!(net, 2, 3)
        add_edge!(net, 3, 4)
        add_edge!(net, 4, 5)

        # In a path graph, all internal vertices are cutpoints
        cuts = cutpoints(net)
        @test Set(cuts) == Set([2, 3, 4])
    end

    # ---------------------------------------------------------------------
    # Missing-dyad policy
    #
    # Every exported measure takes `missing::Symbol=:error`: with masked
    # (unobserved) dyads present it refuses to compute rather than silently
    # reading their face values. `missing=:face` is the explicit opt-in and
    # must reproduce exactly what the same network without a mask returns.
    # ---------------------------------------------------------------------

    # A directed network with two masked dyads: (1,2) is a *present*-face
    # mask (the arc exists in the backing graph) and (4,5) an *absent*-face
    # mask (no arc). `masked` and `observed` have identical structure; only
    # the mask differs.
    function masked_pair(; directed::Bool=true)
        obs = network(6; directed=directed)
        for (i, j) in [(1, 2), (2, 3), (3, 1), (1, 3), (2, 6), (3, 4), (5, 6)]
            add_edge!(obs, i, j)
        end
        masked = copy(obs)
        set_missing_dyad!(masked, 1, 2)   # present face: the arc is there
        set_missing_dyad!(masked, 4, 5)   # absent face: no arc
        return masked, obs
    end

    # (name, thunk) for every exported measure that takes a network. The
    # thunk forwards `missing` so the same call can be made under both
    # policies.
    measure_calls = [
        ("degreecent", (g; kw...) -> degreecent(g; kw...)),
        ("betweenness", (g; kw...) -> betweenness(g; kw...)),
        ("closeness", (g; kw...) -> closeness(g; kw...)),
        ("evcent", (g; kw...) -> evcent(g; kw...)),
        ("bonpow", (g; kw...) -> bonpow(g; exponent=0.1, kw...)),
        ("infocent", (g; kw...) -> infocent(g; kw...)),
        ("flowbet", (g; kw...) -> flowbet(g; kw...)),
        ("centralization(:degree)", (g; kw...) -> centralization(g, :degree; kw...)),
        ("centralization(:betweenness)", (g; kw...) -> centralization(g, :betweenness; kw...)),
        ("centralization(:closeness)", (g; kw...) -> centralization(g, :closeness; kw...)),
        ("centralization(:eigenvector)", (g; kw...) -> centralization(g, :eigenvector; kw...)),
        ("centralization(bonpow)", (g; kw...) -> centralization(g, bonpow; exponent=0.1, kw...)),
        ("gden", (g; kw...) -> gden(g; kw...)),
        ("grecip", (g; kw...) -> grecip(g; kw...)),
        ("grecip(:correlation)", (g; kw...) -> grecip(g; measure=:correlation, kw...)),
        ("transitivity", (g; kw...) -> transitivity(g; kw...)),
        ("transitivity(:local)", (g; kw...) -> transitivity(g; type=:local, kw...)),
        ("gtrans", (g; kw...) -> gtrans(g; kw...)),
        ("gtrans(:strong)", (g; kw...) -> gtrans(g; measure=:strong, kw...)),
        ("dyad_census", (g; kw...) -> dyad_census(g; kw...)),
        ("triad_census", (g; kw...) -> triad_census(g; kw...)),
        ("mutuality", (g; kw...) -> mutuality(g; kw...)),
        ("hierarchy", (g; kw...) -> hierarchy(g; kw...)),
        ("hierarchy(:krackhardt)", (g; kw...) -> hierarchy(g; measure=:krackhardt, kw...)),
        ("efficiency", (g; kw...) -> efficiency(g; kw...)),
        ("connectedness", (g; kw...) -> connectedness(g; kw...)),
        ("component_dist", (g; kw...) -> component_dist(g; kw...)),
        ("reachability", (g; kw...) -> reachability(g; kw...)),
        ("largest_component", (g; kw...) -> largest_component(g; kw...)),
        ("cliques", (g; kw...) -> cliques(g; min_size=2, kw...)),
        ("kcores", (g; kw...) -> kcores(g; kw...)),
        ("cutpoints", (g; kw...) -> cutpoints(g; kw...)),
        ("bicomponents", (g; kw...) -> bicomponents(g; kw...)),
        ("geodist", (g; kw...) -> geodist(g; kw...)),
        ("average_path_length", (g; kw...) -> average_path_length(g; kw...)),
        ("sedist", (g; kw...) -> sedist(g; kw...)),
        ("regular_equivalence", (g; kw...) -> regular_equivalence(g; kw...)),
        ("equiv_clust", (g; kw...) -> equiv_clust(g; k=2, kw...)),
        ("blockmodel", (g; kw...) -> blockmodel(g; k=2, kw...)),
        ("brokerage", (g; kw...) -> brokerage(g, [1, 1, 2, 2, 3, 3]; kw...)),
        ("cug_test", (g; kw...) -> (ct = cug_test(g, gtrans; n_sim=10, rng=MersenneTwister(1), kw...);
                                    (ct.obs_stat, ct.rep_stat))),
        ("gcor", (g; kw...) -> gcor(g, g; kw...)),
        ("gcov", (g; kw...) -> gcov(g, g; kw...)),
        ("qaptest", (g; kw...) -> (qt = qaptest((a, b) -> cor(vec(a), vec(b)),
                                               g, g; n_sim=10,
                                               rng=MersenneTwister(1), kw...);
                                   (qt.testval, qt.dist))),
        ("netlm", (g; kw...) -> netlm(g, g; n_sim=10, rng=MersenneTwister(1), kw...).coefficients),
        # Intercept-only logistic fit: using g as its own predictor completely
        # separates the response, so it has no finite MLE to compare.
        ("netlogit", (g; kw...) -> netlogit(g, AbstractNetwork[]; nullhyp=:classical, kw...).coefficients),
    ]

    @testset "Missing dyads: every measure rejects a masked network" begin
        for directed in (true, false)
            masked, _ = masked_pair(; directed=directed)
            @test n_missing_dyads(masked) == 2
            for (name, f) in measure_calls
                # Default policy is :error: the measure must refuse
                err = try
                    f(masked)
                    nothing
                catch e
                    e
                end
                @test err isa ArgumentError
                # ...and the message must name the routine and the mask
                @test occursin("masked", sprint(showerror, err))
            end
        end
    end

    @testset "Missing dyads: :face reproduces the unmasked answer" begin
        for directed in (true, false)
            masked, obs = masked_pair(; directed=directed)
            for (name, f) in measure_calls
                # The face-value opt-in must return exactly what the same
                # network without a mask returns (the old, documented answer)
                @test isequal(f(masked; missing=:face), f(obs))
            end
        end
    end

    @testset "Missing dyads: present-face and absent-face masks" begin
        masked, obs = masked_pair()

        # The mask does not touch the backing graph: the present-face dyad
        # (1,2) still reads as an edge, the absent-face dyad (4,5) as a non-edge
        @test is_missing_dyad(masked, 1, 2) && has_edge(masked, 1, 2)
        @test is_missing_dyad(masked, 4, 5) && !has_edge(masked, 4, 5)
        @test ne(masked) == ne(obs)

        # Under :face, the present-face masked arc is counted as a tie and the
        # absent-face masked dyad is not — that is precisely the danger the
        # :error default guards against
        @test degreecent(masked; cmode=:outdegree, missing=:face)[1] ==
              degreecent(obs; cmode=:outdegree)[1]
        @test gden(masked; missing=:face) == gden(obs) == 7 / 30
        @test dyad_census(masked; missing=:face) == dyad_census(obs)

        # Masking every dyad of the network still leaves the face values intact
        # but blocks every measure
        allmasked = copy(obs)
        for i in 1:6, j in 1:6
            i == j || set_missing_dyad!(allmasked, i, j)
        end
        @test_throws ArgumentError gden(allmasked)
        @test gden(allmasked; missing=:face) == gden(obs)

        # clear_missing_dyads! declares everything observed again: the default
        # policy works once more
        cleared = clear_missing_dyads!(copy(masked))
        @test n_missing_dyads(cleared) == 0
        @test gden(cleared) == gden(obs)
        @test triad_census(cleared) == triad_census(obs)
    end

    @testset "Missing dyads: policy validation and QAP arguments" begin
        masked, obs = masked_pair()

        # Only :error and :face are policies
        @test Set(MISSING_POLICIES) == Set([:error, :face])
        @test_throws ArgumentError gden(obs; missing=:ignore)
        @test_throws ArgumentError netlm(obs, obs; missing=:ignore, n_sim=10)

        # No SNA measure claims a principled missing-data treatment
        for f in (degree_centrality, density, triad_census, centralization,
                  netlm, netlogit, qaptest)
            @test supports_missing(f) == false
        end

        # A mask on *either* QAP argument is caught: y, a predictor, or the
        # second graph of qaptest
        @test_throws ArgumentError netlm(masked, obs; n_sim=10)
        @test_throws ArgumentError netlm(obs, masked; n_sim=10)
        @test_throws ArgumentError netlogit(obs, [obs, masked]; n_sim=10)
        @test_throws ArgumentError qaptest((a, b) -> cor(vec(a), vec(b)),
                                           obs, masked; n_sim=10)

        # Raw matrices carry no mask, so they are always usable
        A = Float64.(as_matrix(obs))
        fit = netlm(A, A; n_sim=10, rng=MersenneTwister(2))
        @test length(fit.coefficients) == 2

        # The error message names the routine that refused
        msg = sprint(showerror, try
            betweenness(masked)
        catch e
            e
        end)
        @test occursin("betweenness", msg)
        @test occursin("missing=:face", msg)
    end

    @testset "Result metadata protocol (QAP regressions)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)

        # OLS standard errors assume independent dyads; QAP p-values are
        # a separate inference layer. The metadata must preserve that distinction.
        fq = netlm(flo, biz; n_sim=200, rng=Xoshiro(31))
        md = fit_metadata(fq)
        @test md.estimand == :network_regression
        @test md.objective == :least_squares
        @test md.is_exact                      # OLS solves its objective exactly
        @test md.se_method == :ols            # independent-dyad OLS covariance
        # The fit retains the missing-dyad policy it ran under, so the protocol
        # never has to answer ":unspecified" here: an unmasked fit is :none, and
        # a `missing=:face` fit records that it read face values.
        @test md.missing_method == :none
        @test md.tie_method == :not_applicable
        @test any(occursin("QAP permutation p-values", a) for a in md.approximations)
        @test any(occursin("treated as independent observations", a)
                  for a in md.approximations)
        @test any(occursin("standard errors are homoskedastic", a)
                  for a in md.approximations)

        # nullhyp = :classical does no permutation, and says so instead
        fc = netlm(flo, biz; nullhyp=:classical)
        @test !any(occursin("QAP permutation p-values", a)
                   for a in approximations(fc))
        @test any(occursin("nullhyp = :classical", a) for a in approximations(fc))

        # A `missing=:face` fit RECORDS that it read unobserved ties at face
        # value — the opt-in is auditable after the fact, not just at the call
        # site. (This is why the result carries the policy at all.)
        masked = copy(flo)
        set_missing_dyad!(masked, 1, 2)
        @test_throws ArgumentError netlm(masked, biz; nullhyp=:classical)
        ff = netlm(masked, biz; nullhyp=:classical, missing=:face)
        @test missing_method(ff) == :condition_on_face
        @test fit_metadata(ff).missing_method == :condition_on_face
        # ...and a masked PREDICTOR is caught too, not just a masked response
        maskedx = copy(biz)
        set_missing_dyad!(maskedx, 3, 4)
        @test_throws ArgumentError netlm(flo, maskedx; nullhyp=:classical)
        @test missing_method(netlm(flo, maskedx; nullhyp=:classical,
                                   missing=:face)) == :condition_on_face

        # netlogit: an exactly maximized binomial likelihood, with inverse-Fisher
        # standard errors that assume independent dyads — while the p-values come
        # from the permutation null, not from those standard errors
        # Continuous distance covariate supports finite permutation fits; the
        # sparse binary predictor alone can have quasi-separated permutations.
        distance = abs.(FLO_WEALTH .- FLO_WEALTH')
        gl = netlogit(flo, distance; n_sim=200, rng=Xoshiro(32))
        mdl = fit_metadata(gl)
        @test mdl.estimand == :network_logit_regression
        @test mdl.objective == :likelihood
        @test mdl.is_exact
        @test mdl.se_method == :fisher
        @test any(occursin("anticonservative", a) for a in mdl.approximations)
        @test any(occursin("not derived from them", a) for a in mdl.approximations)
    end
end

@testset "Two-mode matrix measures preserve every actor and arc" begin
    for directed in (true, false)
        id = directed ? "bip_dir" : "bip_undir"
        for b in (Network(5; bipartite=2, directed), BipartiteNetwork(2, 3; directed))
            add_edges!(b isa BipartiteNetwork ? b.network : b, [(1,3), (1,4), (2,5)])
            directed && add_edge!(b, 3, 2)
            A = SNA._sociomatrix(b)
            @test size(A) == (5, 5)
            @test A[3,2] == (directed ? 1 : 0)
            full = network_from_matrix(A; directed)
            for (key, value) in (("density", gden(b)),
                                 ("discount_density", gden(b; discount_bipartite=true)),
                                 ("degree", degreecent(b)),
                                 ("bonacich", bonpow(b; exponent=0.05)),
                                 ("flowbet", flowbet(b)))
                @test NetworkCore.check_golden(R_SNA, id * "_" * key, value)
            end
            for f in (betweenness, closeness, evcent, sedist, regular_equivalence)
                @test isequal(f(b), f(full))
            end
            @test isequal(blockmodel(b; k=2).block_matrix, blockmodel(full; k=2).block_matrix)
            # QAP on two-mode data sees only the cross-mode dyads: qaptest's
            # statistic gets the 2 × 3 incidence matrix, and the regression
            # has 2·3 = 6 dyads (plus the 3→2 block when it carries a tie)
            if directed
                @test_throws ArgumentError qaptest((x, y) -> sum(x), b, b; n_sim=5)
                @test netlm(b, full; nullhyp=:classical).n == 12
            else
                shapes = Set{Tuple{Int,Int}}()
                qaptest((x, y) -> (push!(shapes, size(x)); sum(x)), b, b;
                        n_sim=30, threaded=false, rng=Xoshiro(42))
                @test shapes == Set([(2, 3)])
                @test netlm(b, full; nullhyp=:classical).n == 6
            end
        end
    end
end

@testset "Two-mode QAP uses cross-mode dyads only" begin
    # R reference: lm on the 96 cells of an 8 × 12 incidence response and
    # predictor (the cross-mode submatrix), from test/fixtures/r/sna_fuzz.R
    Y = reshape(R_FUZZ.values["twomode_y"], 8, 12)
    X = reshape(R_FUZZ.values["twomode_x"], 8, 12)
    fit = netlm(Y, X; nullhyp=:classical)
    @test fit.n == 96 && !fit.directed
    @test NetworkCore.check_golden(R_FUZZ, "twomode_lm_coefficients", fit.coefficients)
    @test NetworkCore.check_golden(R_FUZZ, "twomode_lm_tstat", fit.tstat)
    @test NetworkCore.check_golden(R_FUZZ, "twomode_gcor", [gcor(Y, X)])
    # The same data as a two-mode Network give the same fit
    B = Network(20; bipartite=8, directed=false)
    for i in 1:8, j in 1:12
        Y[i, j] == 1 && add_edge!(B, i, 8 + j)
    end
    S = zeros(20, 20); S[1:8, 9:20] = X; S[9:20, 1:8] = X'
    @test netlm(B, S; nullhyp=:classical).coefficients ≈ fit.coefficients
    # QAP on rectangular data permutes rows and columns separately
    q = netlm(Y, X; n_sim=50, rng=Xoshiro(3))
    @test size(q.dist) == (50, 2) && q.reps == 50

    # Simulation: independent two-mode networks have a mean slope of about 0
    # (the within-mode zeros used to give 0.25) and the classical test is no
    # longer rejected nine times in ten.
    rng = Xoshiro(20261002)
    slopes = Float64[]
    rejected = 0
    for _ in 1:200
        y = Float64.(rand(rng, 8, 12) .< 0.3)
        x = Float64.(rand(rng, 8, 12) .< 0.3)
        f = netlm(y, x; nullhyp=:classical)
        push!(slopes, f.coefficients[2])
        rejected += f.pgreqabs[2] < 0.05
    end
    @test abs(mean(slopes)) < 3 * std(slopes) / sqrt(200)
    @test rejected / 200 < 0.12
    # A directed two-mode network whose mode-2 → mode-1 block is empty in
    # every argument uses the 1 → 2 block only
    Bd = Network(20; bipartite=8, directed=true)
    for i in 1:8, j in 1:12
        Y[i, j] == 1 && add_edge!(Bd, i, 8 + j)
    end
    Sd = zeros(20, 20); Sd[1:8, 9:20] = X
    @test netlm(Bd, Sd; nullhyp=:classical).n == 96
    @test netlm(Bd, Sd; nullhyp=:classical).coefficients ≈ fit.coefficients
    @test_throws ArgumentError netlm(Y, zeros(9, 12); nullhyp=:classical)
end

@testset "R edge cases and weights" begin
    empty = Network(3)
    for (key, actual) in (("empty_transitivity", transitivity(empty)),
                           ("empty_nonnull_nan", isnan(grecip(empty; measure=:dyadic_nonnull))),
                           ("empty_edgewise_nan", isnan(grecip(empty; measure=:edgewise))),
                           ("singleton_reciprocity_nan", isnan(grecip(Network(1)))))
        @test NetworkCore.check_golden(R_SNA, key, actual)
    end
    @test isnan(grecip(Network(3; directed=false); measure=:edgewise))
    loops = Network(4; loops=true)
    add_edges!(loops, [(1,1), (1,2), (2,3), (3,4)])
    @test NetworkCore.check_golden(R_SNA, "loops_degree", degreecent(loops))
    @test NetworkCore.check_golden(R_SNA, "loops_degree_diag", degreecent(loops; diag=true))
    @test NetworkCore.check_golden(R_SNA, "loops_density", gden(loops))
    @test NetworkCore.check_golden(R_SNA, "loops_density_diag", gden(loops; diag=true))
    W = [0.0 2 1; 1 0 3; 4 2 0]
    weighted = network_from_matrix(W; store_values=true)
    for (key, value) in (("weighted_evcent", evcent(weighted; ignore_eval=false)),
                          ("weighted_bonacich", bonpow(weighted; exponent=0.05, ignore_eval=false)),
                          ("weighted_flowbet", flowbet(weighted; ignore_eval=false)),
                          ("weighted_degree", degreecent(weighted; ignore_eval=false)))
        @test NetworkCore.check_golden(R_SNA, key, value)
    end
    @test_throws ArgumentError degreecent(weighted; cmode=:bogus)
    @test_throws ArgumentError degreecent(weighted; rescale=true, normalized=true)
    @test_throws ArgumentError closeness(weighted; cmode=:bogus)
    @test_throws ArgumentError betweenness(weighted; cmode=:endpoints)
    @test_throws ArgumentError component_dist(weighted; connected=:bogus)
    @test_throws ArgumentError component_dist(weighted; connected=:unilateral)
    @test_throws ArgumentError gtrans(weighted; measure=:rank)
    @test_throws ArgumentError grecip(weighted; measure=:bogus)
    @test_throws ArgumentError sedist(weighted; method=:bogus)
    @test_throws ArgumentError cutpoints(weighted; connected=:bogus)
    @test_throws ArgumentError cliques(weighted; symmetrize=:bogus)
    # Bipartite adjacency has equally large ±rho; direct eigen solve must obey
    # the eigenvector equation, unlike unshifted power iteration on this path.
    path = Network(3; directed=false); add_edges!(path, [(1,2),(2,3)])
    ev = evcent(path)
    @test SNA._sociomatrix(path) * ev ≈ sqrt(2) * ev
    @test isempty(largest_component(Network(0)))
end

@testset "StatsAPI and scheduling-independent QAP" begin
    flo, biz = florentine(), load_dataset(:florentine_business)
    d = abs.(FLO_WEALTH .- FLO_WEALTH')
    q1 = qaptest(gcor, flo, biz; n_sim=30, rng=Xoshiro(2), threaded=false)
    q2 = qaptest(gcor, flo, biz; n_sim=30, rng=Xoshiro(2), threaded=true)
    @test q1.dist == q2.dist
    for fitfun in (netlm, netlogit), nullhyp in (:qapy, :qapx, :qapspp)
        # Continuous predictors keep this scheduling/StatsAPI test identified;
        # separated binary permutations have a separate refusal regression below.
        predictors = fitfun === netlogit ? [d, FLO_WEALTH .+ FLO_WEALTH'] : [biz, d]
        a = fitfun(flo, predictors; n_sim=30, rng=Xoshiro(44), threaded=false, nullhyp)
        b = fitfun(flo, predictors; n_sim=30, rng=Xoshiro(44), threaded=true, nullhyp)
        @test a.dist == b.dist
        @test a.pgreqabs == b.pgreqabs
        @test coef(a) === a.coefficients
        @test stderror(a) ≈ sqrt.(diag(vcov(a)))
        @test size(confint(a)) == (3,2)
        @test size(vcov(a)) == (3,3)
        @test isposdef(Symmetric(vcov(a)))
        @test nobs(a) == 120
        @test aic(a) ≈ -2loglikelihood(a) + 2dof(a)
        @test bic(a) ≈ -2loglikelihood(a) + log(nobs(a))*dof(a)
        @test coeftable(a) isa NetworkCore.CoefficientTable
        @test all(values(NetworkCore.check_statsapi(a; strict=true)))
        # `coefnames` is R's names(coef(fit)), the table's names, and a copy
        @test coefnames(a) == coeftable(a).names == a.names
        @test length(coefnames(a)) == 3 && coefnames(a)[1] == "(intercept)"
        @test coefnames(a) !== a.names
        @test SNA.coefnames === NetworkCore.coefnames === NetworkCore.StatsAPI.coefnames
        @test all(values(NetworkCore.check_statsapi(a;
            required=(NetworkCore.STATSAPI_VERBS..., :coefnames), strict=true)))
    end
    zero_tail = QAPTestResult(2.0, [0.0,1.0], 0.0, 1.0, 2)
    @test occursin("<0.5", sprint(show, zero_tail))
    # Imported fixture/optimizer helpers stay qualified in the public namespace.
    @test :Network in names(SNA)
    @test !(:load_golden in names(SNA))
    @test !(:bootstrap_cov in names(SNA))
    @test !(:NetworkCore in names(SNA))
    @test Network(5) isa NetworkCore.Network
    for f in (gden, degreecent, betweenness, netlm, netlogit, component_dist, gcor)
        @test missing_policies(f) == (:error, :face)
    end
end


@testset "Regular equivalence is invariant to actor ordering" begin
    net = Network(6)
    add_edges!(net, [(1,2), (1,3), (2,4), (3,5), (5,6), (6,5)])
    A = SNA._sociomatrix(net)
    p = [6,3,1,5,2,4]
    @test regular_equivalence(network_from_matrix(A[p,p])) ≈ regular_equivalence(net)[p,p]
    @test_throws ArgumentError regular_equivalence(net; ignore_eval=false)
    @test_throws ArgumentError regular_equivalence(net; maxiter=0)
    # the development-time keyword `max_iter` is `maxiter`; there is no alias
    @test_throws MethodError regular_equivalence(net; max_iter=50)
    flo = florentine()
    @test_throws ArgumentError netlm(flo, [flo, flo]; nullhyp=:classical)
    # y regressed on itself is completely separated: warned, not refused
    selfit = @test_logs (:warn, r"separation") netlogit(flo, flo; nullhyp=:classical)
    @test selfit.separated == ["x1"] && !selfit.converged
end


@testset "SNA alone exports its constructor" begin
    code = "using SNA; @assert Network(5) isa Network; @assert !isdefined(Main, :load_golden)"
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(pkgdir(SNA)) -e $code`
    @test success(pipeline(cmd; stdout=devnull))
end


@testset "Quasi-complete logistic separation" begin
    response = zeros(4, 4)
    predictor = zeros(4, 4)
    for edge in ((1, 2), (2, 1), (3, 4), (4, 3))
        response[edge...] = 1
    end
    predictor[1, 2] = predictor[2, 1] = 1
    # Predictor=1 has only successes; predictor=0 contains both outcomes.
    # Its slope tends to +Inf even though an optimizer can stop numerically.
    # The ecosystem's separation policy: warn, converged == false, the
    # separated term flagged, inference withheld (NaN), for every null
    # hypothesis and at any scale of the predictor.
    for (yy, xx, nullhyp) in ((response, predictor, :classical),
                              (response, 1e6 .* predictor, :classical),
                              (1 .- response - Matrix{Float64}(I, 4, 4), predictor, :classical),
                              (response, predictor, :qapspp),
                              (response, predictor, :qapy))
        fit = @test_logs (:warn, r"netlogit: the maximum-likelihood estimate does not exist \(separation\)") netlogit(
            yy, xx; nullhyp, n_sim=20, rng=Xoshiro(1))
        @test !fit.converged && !is_exact(fit)
        @test fit.separated == ["x1"] && fit.separation.separated && fit.separation.certified
        @test all(isnan, fit.tstat) && all(isnan, fit.pgreqabs) && all(isnan, fit.pleeq)
        @test all(isnan, confint(fit))
        @test fit.dist === nothing && fit.reps == 0
        @test all(isnan, coeftable(fit).p_values)
        @test any(occursin("separation: the coefficient(s) on `x1`", a) for a in approximations(fit))
        @test occursin("Warning: separation", sprint(show, fit))
        @test isfinite(coef(fit)[1])                 # the estimates stay, for diagnosis
    end
    # Adding a failure to the predictor=1 stratum identifies a finite slope.
    predictor[1, 3] = 1
    fitted = @test_logs netlogit(response, predictor; nullhyp=:classical)
    @test fitted.converged && isempty(fitted.separated) && !fitted.separation.separated
    @test coef(fitted)[2] ≈ log(2 / 1) - log(2 / 7) atol=1e-7
end


@testset "Separation verdict on exact margins; Wald replicates withheld, not refused" begin
    # The verdict is NetworkCore's (an exact certificate on the Float64
    # data). A tiny interior reversal still gives overlap in exact
    # arithmetic, but the overlap is at the 1e-16 level: the design is
    # separated up to the rounding of its entries, and the verdict says so
    # without claiming a certificate.
    y = [0.0, 1.0, 0.0, 1.0]
    X = hcat(ones(4), [-1.0, -1e-16, 1e-16, 1.0])
    v = NetworkCore.logistic_separation(X, Bool.(y))
    @test v.separated && !v.certified
    w = NetworkCore.logistic_separation(X, Bool[0, 0, 1, 1])
    @test w.separated && w.certified && w.terms == [2]
    # Overlap at x = 0 alone does not stop the slope: the x = ±1 rows are
    # still predicted perfectly (quasi-complete separation)
    Q = hcat(ones(4), [-1.0, 0.0, 0.0, 1.0])
    vq = NetworkCore.logistic_separation(Q, Bool.(y))
    @test vq.separated && vq.certified && vq.terms == [2] && vq.units == [1, 4]
    # overlap at both ends: a finite maximum
    O = hcat(ones(4), [-1.0, -1.0, 1.0, 1.0])
    @test !NetworkCore.logistic_separation(O, Bool.(y)).separated
    @test SNA._logit_fit(O, y).converged
    flo, biz = florentine(), load_dataset(:florentine_business)
    # sna's Wald statistic does not exist on a separated permutation: with
    # statistic=:wald the permutation p-value of the business coefficient is
    # withheld (NaN, warned), identically in serial and threaded runs, and the
    # intercept's, whose replicates are all defined, is not
    fits = map((false, true)) do threaded
        @test_logs (:warn, r"permutation replicates are separated") netlogit(
            flo, biz; statistic=:wald, n_sim=200, rng=Xoshiro(32), threaded)
    end
    @test isequal(fits[1].dist, fits[2].dist)
    for f in fits
        @test f.converged && isempty(f.separated)         # the observed design is fine
        @test isnan(f.pgreqabs[2]) && isnan(f.pleeq[2]) && isnan(f.pgreq[2])
        @test count(isnan, f.dist[:, 2]) > 0
        @test any(occursin("withheld", a) for a in approximations(f))
    end
end

@testset "netlogit QAP runs on separated permutations" begin
    # The canonical marriage ~ business QAP: about 6 % of permutations are
    # separated. The default signed-root LR statistic exists on all of them,
    # so every replicate counts; sna reports p = 0.000 / 0.001.
    flo, biz = florentine(), load_dataset(:florentine_business)
    for nullhyp in (:qapspp, :qapy, :qapx)
        fits = [netlogit(flo, biz; nullhyp, n_sim=300, rng=Xoshiro(32), threaded)
                for threaded in (false, true)]
        @test fits[1].dist == fits[2].dist            # thread-count independent
        f = fits[1]
        @test all(isfinite, f.dist)
        @test size(f.dist) == (300, 2)
        @test f.pgreqabs[2] <= 0.01                   # business ties
    end
    # At least one of those replicates is a separated design (the reason the
    # Wald statistic fails): the permuted business ties of seed 32's 22nd draw
    Y = SNA._sociomatrix(flo); B = SNA._sociomatrix(biz)
    idx = SNA._dyad_indices(16, false)
    yv = SNA._gvectorize(Y, idx)
    seps = count(1:200) do s
        X = hcat(ones(120), SNA._gvectorize(SNA._rmperm(Xoshiro(s), B), idx))
        NetworkCore.logistic_separation(X, Bool.(yv)).separated
    end
    @test seps > 0
    # Size under the null: y independent of x, DSP QAP with the LR statistic
    rng = Xoshiro(99)
    pvals = Float64[]
    for _ in 1:60
        n = 10
        y = Float64.(rand(rng, n, n) .< 0.25)
        x = Float64.(rand(rng, n, n) .< 0.15)
        for i in 1:n
            y[i, i] = x[i, i] = 0
        end
        f = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            netlogit(y, x; n_sim=99, rng=Xoshiro(rand(rng, UInt)), mode=:digraph)
        end
        isempty(f.separated) || continue            # observed design separated: no p-value
        push!(pvals, f.pgreqabs[2])
    end
    @test length(pvals) >= 50
    @test mean(pvals .<= 0.05) <= 0.15
    @test 0.3 <= mean(pvals) <= 0.7                  # roughly uniform
end


@testset "Binary-stratum separation agrees with analytical likelihood" begin
    for n0 in (2, 3), n1 in (2, 3), successes0 in 0:n0, successes1 in 0:n1
        X = hcat(ones(n0+n1), vcat(zeros(n0), ones(n1)))
        y = vcat(ones(successes0), zeros(n0-successes0),
                 ones(successes1), zeros(n1-successes1))
        if successes0 in (0,n0) || successes1 in (0,n1)
            fit = SNA._logit_fit(X, y)
            @test fit.verdict.separated && fit.verdict.certified && !fit.converged
            # A degenerate x = 1 stratum is carried by the slope alone; a
            # degenerate x = 0 stratum needs the intercept, with the slope
            # moving opposite so that the x = 1 stratum is unchanged.
            deg0 = successes0 in (0, n0); deg1 = successes1 in (0, n1)
            deg1 && !deg0 && @test fit.verdict.terms == [2]
            deg0 && !deg1 && @test fit.verdict.terms == [1, 2]
            deg0 && deg1 && @test !isempty(fit.verdict.terms)
        else
            fit = SNA._logit_fit(X, y)
            alpha = log(successes0 / (n0-successes0))
            beta = log(successes1 / (n1-successes1)) - alpha
            @test fit.coef ≈ [alpha, beta] atol=1e-7
            @test fit.converged
        end
    end
end


# ---------------------------------------------------------------------------
# sna 2.8 fuzz corpus (test/fixtures/r/sna_fuzz.R): 60 seeded random networks,
# directed and undirected, with self-loops, isolates and disconnected parts.
# Every value is R's; NaN must match NaN (sna's undefined ratios).
# ---------------------------------------------------------------------------
@testset "Fuzz corpus: every measure agrees with sna 2.8 (and igraph), Int and Int32 ids" begin
    same(e, a, atol) = length(e) == length(a) &&
        all((isnan(x) && isnan(y)) || isapprox(x, y; atol) for (x, y) in zip(e, a))
    vals = R_FUZZ.values
    ids = sort(unique(parse(Int, m.captures[1])
                      for k in keys(vals) for m in (match(r"^g(\d+)_meta$", k),) if m !== nothing))
    @test length(ids) == 60
    # The whole corpus runs twice: with the default Int vertex ids and with
    # Int32 ids, which NetworkCore supports. Six measures used to throw a
    # MethodError on Int32 through `n::Int` signatures of internal helpers.
    for T in (Int, Int32)
    compared = Dict{String,Int}()
    for g in ids
        n, directed, loops = Int.(vals["g$(g)_meta"])
        net = Network{T}(; n=n, directed=directed == 1, loops=loops == 1)
        E = Int.(vals["g$(g)_edges"])
        for k in 1:2:length(E)
            add_edge!(net, E[k], E[k+1])
        end
        calls = [
            "degree" => () -> degreecent(net),
            "degree_in" => () -> degreecent(net; cmode=:indegree),
            "degree_out" => () -> degreecent(net; cmode=:outdegree),
            "betweenness" => () -> betweenness(net),
            "betweenness_undirected" => () -> betweenness(net; cmode=:undirected),
            "closeness" => () -> closeness(net),
            "closeness_suminv" => () -> closeness(net; cmode=:suminvdir),
            "closeness_gil_schmidt" => () -> closeness(net; cmode=:gil_schmidt),
            "bonpow" => () -> bonpow(net; exponent=0.3),
            "infocent" => () -> infocent(net),
            "flowbet" => () -> flowbet(net),
            "centralization_degree" => () -> centralization(net, degreecent),
            "centralization_betweenness" => () -> centralization(net, betweenness),
            "centralization_closeness" => () -> centralization(net, closeness),
            "gden" => () -> gden(net),
            "grecip_dyadic" => () -> grecip(net),
            "grecip_dyadic_nonnull" => () -> grecip(net; measure=:dyadic_nonnull),
            "grecip_edgewise" => () -> grecip(net; measure=:edgewise),
            "grecip_edgewise_lrr" => () -> grecip(net; measure=:edgewise_lrr),
            "grecip_correlation" => () -> grecip(net; measure=:correlation),
            "gtrans_weak" => () -> gtrans(net),
            "gtrans_strong" => () -> gtrans(net; measure=:strong),
            "gtrans_weakcensus" => () -> gtrans(net; measure=:weakcensus),
            "gtrans_strongcensus" => () -> gtrans(net; measure=:strongcensus),
            "gtrans_correlation" => () -> gtrans(net; measure=:correlation),
            "dyad_census" => () -> collect(values(dyad_census(net))),
            "triad_census" => () -> triad_census(net),
            "hierarchy" => () -> hierarchy(net),
            "hierarchy_krackhardt" => () -> hierarchy(net; measure=:krackhardt),
            "efficiency" => () -> efficiency(net),
            "connectedness" => () -> connectedness(net),
            "geodist" => () -> vec(geodist(net).gdist),
            "geodist_counts" => () -> vec(geodist(net).counts),
            "kcores_freeman" => () -> kcores(net),
            "kcores_indegree" => () -> kcores(net; cmode=:indegree),
            "kcores_outdegree" => () -> kcores(net; cmode=:outdegree),
            "local_clustering_weak" => () -> transitivity(net; type=:local, cmode=:weak),
            "average_clustering_weak" => () -> transitivity(net; type=:average, cmode=:weak),
            "gcor_transpose" => () -> gcor(net, Matrix(transpose(SNA._sociomatrix(net)));
                                           mode=is_directed(net) ? :digraph : :graph),
            "gcov_transpose" => () -> gcov(net, Matrix(transpose(SNA._sociomatrix(net)));
                                           mode=is_directed(net) ? :digraph : :graph),
        ]
        for cn in (:strong, :weak, :recursive), part in (:membership, :csize, :cdist)
            push!(calls, "components_$(cn)_$(part)" => () -> getfield(component_dist(net; connected=cn), part))
        end
        for m in (:hamming, :correlation, :euclidean, :gamma, :exact)
            push!(calls, "sedist_$m" => () -> vec(sedist(net; method=m)))
        end
        for (key, f) in calls
            haskey(vals, "g$(g)_$key") || continue        # sna refused this one
            actual = Float64.(vec(collect(f())))
            ok = same(Float64.(vals["g$(g)_$key"]), actual, 1e-9)
            ok || @info "fuzz mismatch" T g key actual vals["g$(g)_$key"]
            @test ok
            compared[key] = get(compared, key, 0) + 1
        end
        if T === Int32
            # the two path-based routines the fixture does not record: the
            # Int32 network gives exactly the Int network's answer
            net64 = network(n; directed=directed == 1, loops=loops == 1)
            for k in 1:2:length(E)
                add_edge!(net64, E[k], E[k+1])
            end
            @test isequal(average_path_length(net), average_path_length(net64))
            @test layout_kamada_kawai(net) == layout_kamada_kawai(net64)
            # cug_test draws its null graphs from the size, not the id type
            for cmode in (:size, :edges, :dyad_census)
                c32 = cug_test(net, gden; cmode, n_sim=3, rng=Xoshiro(g), threaded=false)
                c64 = cug_test(net64, gden; cmode, n_sim=3, rng=Xoshiro(g), threaded=false)
                @test c32.obs_stat == c64.obs_stat && c32.rep_stat == c64.rep_stat
            end
        end
    end
    @test length(compared) == 54      # every measure key of the fixture was compared
    @test all(>=(40), values(compared))
    end
end

@testset "Measure edge cases follow sna" begin
    # Bonacich power on a well-conditioned system whose determinant underflows:
    # a perfect matching of 200 actors, β = 0.5,
    # condition number 3, det(I − βA) = 3e-13. sna gives all ones.
    m = network(200; directed=false)
    for i in 1:2:199
        add_edge!(m, i, i + 1)
    end
    @test NetworkCore.check_golden(R_FUZZ, "matching_bonpow", bonpow(m; exponent=0.5))
    # A singular system is still refused (NaN with a warning)
    star = network(3; directed=false); add_edges!(star, [(1, 2), (1, 3)])
    @test all(isnan, @test_logs (:warn, r"singular") bonpow(star; exponent=1 / sqrt(2)))
    # No ties: no scale α, as sna
    @test all(isnan, bonpow(network(3)))

    # evcent of a nilpotent (acyclic) digraph has no Perron vector: NaN
    dag = network(6)
    add_edges!(dag, [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (1, 6)])
    @test all(isnan, evcent(dag))
    @test all(isnan, evcent(network(4)))
    ring = network(3); add_edges!(ring, [(1, 2), (2, 3), (3, 1)])
    @test evcent(ring) ≈ fill(1 / sqrt(3), 3)

    # Self-loops: efficiency ignores them, kcores too
    looped = network(3; loops=true)
    add_edges!(looped, [(1, 2), (2, 3), (1, 1), (2, 2)])
    plain = network(3); add_edges!(plain, [(1, 2), (2, 3)])
    @test efficiency(looped) == efficiency(plain) == 1.0
    @test kcores(looped) == kcores(plain)
    @test efficiency(looped) >= 0
    # Undefined ratios are NaN, as in sna
    @test isnan(efficiency(network(3)))
    @test isnan(hierarchy(network(1)))
    @test isnan(hierarchy(network(3); measure=:krackhardt))
    @test connectedness(network(1)) == 1.0
    # The development-time names were never released and have no aliases
    for old in (:bonacich_power, :geodesic_distance, :reciprocity,
                :structural_equivalence, :consensus, :components)
        @test !isdefined(SNA, old)
    end
end

@testset "consensus_clustering is order-invariant and transitive" begin
    parts = [[1, 1, 2, 2, 3], [1, 1, 1, 2, 3], [2, 2, 1, 1, 3]]
    labels = consensus_clustering(parts)
    @test labels == [1, 1, 2, 2, 3]
    p = [5, 3, 1, 4, 2]
    permuted = consensus_clustering([c[p] for c in parts])
    # same partition, up to relabelling
    @test [permuted[i] == permuted[j] for i in 1:5, j in 1:5] ==
          [labels[p][i] == labels[p][j] for i in 1:5, j in 1:5]
    # a chain of pairwise agreements is one cluster (transitive)
    @test consensus_clustering([[1, 1, 2], [2, 1, 1]]; threshold=0.5) == [1, 1, 1]
    @test_throws ArgumentError consensus_clustering(Vector{Int}[])
    @test_throws ArgumentError consensus_clustering([[1, 2], [1, 2, 3]])
end

@testset "layout_kamada_kawai minimizes the Kamada–Kawai stress" begin
    net = network(8; directed=false)
    add_edges!(net, [(1, 2), (2, 3), (3, 4), (4, 1), (4, 5), (5, 6), (6, 7), (7, 8)])
    D = geodist(net).gdist
    stress(X) = sum((norm(X[i, :] - X[j, :]) - D[i, j])^2 / D[i, j]^2
                    for i in 1:8 for j in (i+1):8)
    kk = layout_kamada_kawai(net)
    mds = SNA._classical_mds(D)
    @test stress(kk) < stress(mds)                 # better than its MDS start
    @test kk == layout_kamada_kawai(net)           # deterministic
    # a path is laid out on a line at unit spacing
    path = network(4; directed=false); add_edges!(path, [(1, 2), (2, 3), (3, 4)])
    P = layout_kamada_kawai(path)
    @test [norm(P[i, :] - P[i+1, :]) for i in 1:3] ≈ ones(3) atol = 1e-4
end

@testset "Keyword vocabulary (ecosystem convention)" begin
    # `maxiter`, `n_sim`, `rng`; never `max_iter`, `reps`, `seed`
    for f in (qaptest, netlm, netlogit, regular_equivalence, layout_kamada_kawai,
              degreecent, betweenness, closeness, evcent, bonpow, rgraph, rgnp, rgnm)
        kws = Set{Symbol}()
        for m in methods(f)
            union!(kws, Base.kwarg_decl(m))
        end
        @test :seed ∉ kws
        @test :n_sims ∉ kws
        @test :reps ∉ kws && :max_iter ∉ kws      # no development-time aliases
    end
    @test :n_sim in Base.kwarg_decl(first(methods(qaptest)))
    @test :maxiter in Base.kwarg_decl(first(methods(regular_equivalence)))
end

@testset "Every export is documented with a runnable example" begin
    # Execute every ```julia block of every docstring SNA itself owns, in a
    # fresh module (examples must be self-contained programs). Re-exported
    # NetworkCore/StatsAPI names are documented (and tested) where they live.
    own = Dict{Symbol,String}()
    for (binding, multidoc) in Base.Docs.meta(SNA)
        own[binding.var] = join((join(string.(d.text), "\n") for d in values(multidoc.docs)), "\n")
    end
    for n in names(SNA)
        (n === :SNA || !isdefined(SNA, n)) && continue
        Base.which(SNA, n) === SNA || continue       # re-exports are documented upstream
        @test haskey(own, n)
        haskey(own, n) || continue
        @test occursin("# Example", own[n])
        blocks = [m.captures[1] for m in eachmatch(r"```julia\n(.*?)```"s, own[n])]
        @test !isempty(blocks)
        for code in blocks
            mod = Module(gensym(string(n)))
            ok = try
                Base.include_string(mod, code)
                true
            catch e
                @error "docstring example of $n failed" exception = (e, catch_backtrace())
                false
            end
            @test ok
        end
    end
end

@testset "Julia hygiene: Aqua and method ambiguities" begin
    @test isempty(Test.detect_ambiguities(SNA; recursive=true))
    Aqua.test_all(SNA)
end


# ---------------------------------------------------------------------------
# sna 2.8 references for brokerage, equiv.clust/cutree/blockmodel and
# cug.test (test/fixtures/r/sna_inference.R)
# ---------------------------------------------------------------------------
function inference_graph(g)
    n, directed = Int.(R_INF.values["g$(g)_meta"])
    net = network(n; directed=directed == 1)
    E = Int.(R_INF.values["g$(g)_edges"])
    for k in 1:2:length(E)
        add_edge!(net, E[k], E[k+1])
    end
    return net
end
same_values(e, a; atol=1e-9) = length(e) == length(a) &&
    all((isnan(x) && isnan(y)) || x == y || isapprox(x, y; atol) for (x, y) in zip(e, a))

@testset "brokerage agrees with sna::brokerage (Gould–Fernandez)" begin
    for g in 1:24
        net = inference_graph(g)
        cl = Int.(R_INF.values["g$(g)_cl"])
        b = brokerage(net, cl)
        for key in (:raw_nli, :exp_nli, :sd_nli, :z_nli, :raw_gli, :exp_gli, :sd_gli,
                    :z_gli, :exp_grp, :sd_grp, :clid, :n)
            @test same_values(Float64.(R_INF.values["g$(g)_brokerage_$key"]),
                              Float64.(vec(collect(getfield(b, key)))))
        end
    end
    # Role definitions on i → j → k, by the classes of (i, j, k)
    chain = network(3); add_edges!(chain, [(1, 2), (2, 3)])
    for (cl, role) in (([1, 1, 1], 1), ([1, 2, 1], 2), ([1, 1, 2], 3), ([1, 2, 2], 4), ([1, 2, 3], 5))
        @test brokerage(chain, cl).raw_nli[2, :] == [k == role || k == 6 ? 1.0 : 0.0 for k in 1:6]
    end
    add_edge!(chain, 1, 3)                      # a direct tie: nothing to broker
    @test all(iszero, brokerage(chain, [1, 1, 1]).raw_nli)
    # classes from a vertex attribute; loops ignored; wrong length refused
    net = inference_graph(1)
    cl = Int.(R_INF.values["g1_cl"])
    for v in 1:nv(net)
        set_vertex_attribute!(net, :grp, v, cl[v])
    end
    @test brokerage(net, :grp).raw_nli == brokerage(net, cl).raw_nli
    @test brokerage(net, :grp).roles == ["w_I", "w_O", "b_IO", "b_OI", "b_O", "t"]
    @test_throws ArgumentError brokerage(net, cl[1:2])
end

@testset "equiv_clust defaults are sna::equiv.clust + cutree" begin
    for g in 1:24
        net = inference_graph(g)
        for k in 2:4
            @test equiv_clust(net; k) == Int.(R_INF.values["g$(g)_equiv_clust_k$k"])
        end
        for cm in (:average, :single)
            @test equiv_clust(net; k=3, method=:euclidean, cluster_method=cm) ==
                  Int.(R_INF.values["g$(g)_equiv_clust_euclidean_$(cm)_k3"])
        end
        @test same_values(R_INF.values["g$(g)_blockmodel_k3"], vec(blockmodel(net; k=3).block_matrix))
    end
    net = inference_graph(2)
    # The development-time spellings `method=:structural`/`:regular` were
    # never released; they are refused as any unknown sedist method is
    @test_throws ArgumentError equiv_clust(net; k=3, method=:structural)
    @test_throws ArgumentError equiv_clust(net; k=3, method=:regular)
    @test_throws ArgumentError equiv_clust(net; cluster_method=:ward)
    @test_throws ArgumentError equiv_clust(net; equiv_fun=gden)
    @test_throws ArgumentError equiv_clust(net; method=:bogus)
end

@testset "cug_test agrees with sna::cug.test" begin
    flo, samp = florentine(), sampson_like()
    J = 4000
    R = R_INF.values["cug_reps"][1]
    cases = [("flo_gtrans", flo, g -> gtrans(g)), ("samp_gtrans", samp, g -> gtrans(g)),
             ("samp_grecip", samp, g -> grecip(g; measure=:edgewise))]
    # Seeds are fixed integers: `hash` of a string differs across Julia
    # versions, which would give each CI cell a different stream
    for (ci, (name, net, f)) in enumerate(cases), (mi, cmode) in enumerate((:size, :edges, :dyad_census))
        key = "cug_$(name)_$(cmode)"
        ct = cug_test(net, f; cmode, n_sim=J, rng=Xoshiro(1000 + 10ci + mi))
        @test ct isa CUGTestResult && ct.reps == J && ct.cmode == cmode
        @test ct.obs_stat ≈ R_INF.values[key * "_obs"][1] atol = 1e-12
        # Null distribution within 5 Monte Carlo standard errors of R's
        sd = R_INF.values[key * "_sd"][1]
        se = sqrt(1 / R + 1 / J)
        @test abs(mean(ct.rep_stat) - R_INF.values[key * "_mean"][1]) <= 5 * sd * se + 1e-12
        @test abs(std(ct.rep_stat) - sd) <= 0.1 * sd + 1e-12
        for (tail, val) in (("pgteobs", ct.pgteobs), ("plteobs", ct.plteobs))
            p = R_INF.values[key * "_" * tail][1]
            @test abs(val - p) <= 5 * sqrt(p * (1 - p)) * se + 2 / J
        end
    end
    # The conditioning holds exactly in every replicate
    dc = dyad_census(samp)
    e = cug_test(samp, g -> ne(g); cmode=:edges, n_sim=50, rng=Xoshiro(1))
    @test all(==(ne(samp)), e.rep_stat) && e.pgteobs == e.plteobs == 1.0
    for field in (:mutual, :asymmetric, :null)
        d = cug_test(samp, g -> getfield(dyad_census(g), field); cmode=:dyad_census,
                     n_sim=50, rng=Xoshiro(2))
        @test all(==(getfield(dc, field)), d.rep_stat)
    end
    u = cug_test(flo, g -> is_directed(g) ? -1 : ne(g); cmode=:dyad_census, n_sim=20, rng=Xoshiro(3))
    @test all(==(ne(flo)), u.rep_stat)             # undirected stays undirected
    # :size draws ties with probability 1/2
    sz = cug_test(samp, gden; cmode=:size, n_sim=400, rng=Xoshiro(4))
    @test abs(mean(sz.rep_stat) - 0.5) < 0.01
    # Uniformity over dyad-census-conditioned digraphs: n = 3, M = 1, A = 1 has
    # 3 · 2 · 2 = 12 equally likely digraphs
    tiny = network(3); add_edges!(tiny, [(1, 2), (2, 1), (2, 3)])
    code(g) = sum(2.0^(3 * (i - 1) + (j - 1)) for i in 1:3, j in 1:3 if has_edge(g, i, j))
    draws = cug_test(tiny, code; cmode=:dyad_census, n_sim=6000, rng=Xoshiro(5)).rep_stat
    freq = [count(==(c), draws) for c in unique(draws)]
    @test length(freq) == 12
    @test all(f -> abs(f - 500) < 5 * sqrt(500 * 11 / 12), freq)
    # Reproducible, thread-count independent, keyword passing, loops ignored
    a = cug_test(samp, gtrans; cmode=:edges, n_sim=40, rng=Xoshiro(9), threaded=false)
    b = cug_test(samp, gtrans; cmode=:edges, n_sim=40, rng=Xoshiro(9), threaded=true)
    @test a.rep_stat == b.rep_stat
    @test cug_test(samp, grecip; measure=:edgewise, n_sim=10, rng=Xoshiro(1)).obs_stat ==
          grecip(samp; measure=:edgewise)
    looped = network(4; loops=true); add_edges!(looped, [(1, 1), (1, 2), (2, 3)])
    @test cug_test(looped, g -> ne(g); cmode=:edges, n_sim=5, rng=Xoshiro(1)).obs_stat == 2
    @test occursin("Conditional Uniform Graph", sprint(show, a))
    @test_throws ArgumentError cug_test(samp, gtrans; cmode=:bogus)
    @test_throws ArgumentError cug_test(samp, gtrans; n_sim=0)
    @test_throws ArgumentError cug_test(Network(5; bipartite=2), gden)
end
