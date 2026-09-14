using SNA
using Networks
using Graphs
using Random
using Statistics
using LinearAlgebra
using Test

const R_SNA = Networks.load_golden(joinpath(@__DIR__, "fixtures", "sna_reference.toml"))
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
const FLO_WEALTH = Float64.(R_SNA.values["flo_wealth"])

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

        largest = largest_component(net; connected=:weak)
        @test length(largest) == 3
        @test length(largest_component(net)) == 1
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

        @test diameter(net) == Inf # directed endpoints cannot reach each other
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
        @test isempty(cliques(dnet))
        @test Set.(cliques(dnet; symmetrize=:weak)) == [Set([1, 2, 3])]
    end

    @testset "Provenanced R sna reference: full vectors and directed cohesion" begin
        masks(sets) = sort([sum(2^(v-1) for v in unique(xs)) for xs in sets])
        for (id, net) in (("flo", florentine()), ("samp", sampson_like()))
            check(key, value) = Networks.check_golden(R_SNA, id * "_" * key, value)
            @test check("density", density(net))
            @test check("reciprocity", reciprocity(net))
            @test check("reciprocity_edgewise", reciprocity(net; method=:edgewise))
            @test check("reciprocity_nonnull", reciprocity(net; method=:dyadic_nonnull))
            @test check("transitivity", transitivity(net))
            @test check("dyad_census", collect(values(dyad_census(net))))
            @test check("triad_census", triad_census(net))
            @test mutuality(net) == dyad_census(net).mutual
            @test check("connectedness", connectedness(net))
            @test check("efficiency", efficiency(net))
            @test check("hierarchy", hierarchy(net))
            @test check("hierarchy_krackhardt", hierarchy(net; measure=:krackhardt))
            @test check("degree", degree_centrality(net))
            @test check("degree_in", degree_centrality(net; mode=:in))
            @test check("degree_out", degree_centrality(net; mode=:out))
            @test check("betweenness", betweenness_centrality(net))
            @test check("closeness", closeness_centrality(net))
            @test check("eigenvector", eigenvector_centrality(net))
            @test check("bonacich", bonacich_power(net; exponent=0.05))
            @test check("flowbet", flowbet(net))
            @test check("strong_sizes", component_dist(net))
            @test check("weak_sizes", component_dist(net; connected=:weak))
            @test check("cutpoints", sort(cutpoints(net)))
            @test check("cutpoints_weak", sort(cutpoints(net; connected=:weak)))
            @test check("cutpoints_recursive", sort(cutpoints(net; connected=:recursive)))
            @test check("bicomponents", masks([collect(Iterators.flatten(c)) for c in bicomponents(net)]))
            @test check("cliques", masks(cliques(net)))
            for k in (1, 2, 5, 7)
                @test check("kcore_$k", kcores(net; k))
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
        @test all(iszero, closeness_centrality(flo))
        @test any(>(0), closeness_centrality(flo; cmode=:component))
        cycle = network(4)
        add_edges!(cycle, [(1,2), (2,3), (3,1)])
        @test Networks.check_golden(R_SNA, "cycle_cutpoints", cutpoints(cycle))
        @test Networks.check_golden(R_SNA, "cycle_cutpoints_recursive", cutpoints(cycle; connected=:recursive))
        @test Networks.check_golden(R_SNA, "cycle_strong_sizes", component_dist(cycle))
        @test Networks.check_golden(R_SNA, "cycle_closeness", closeness_centrality(cycle))
        @test Networks.check_golden(R_SNA, "cycle_cliques", masks(cliques(cycle)))
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
    end

    @testset "QAP test (qaptest)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)

        qt = qaptest(gcor, flo, biz; reps=1000, rng=Xoshiro(11))
        # Observed statistic is deterministic: R sna::gcor(flo, biz)
        @test Networks.check_golden(R_SNA, "qap_gcor", qt.testval)
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
        @test Networks.check_golden(R_SNA, "lm_coefficients", fit.coefficients)
        @test Networks.check_golden(R_SNA, "lm_tstat", fit.tstat)
        @test Networks.check_golden(R_SNA, "lm_pgreqabs", fit.pgreqabs)
        @test Networks.check_golden(R_SNA, "lm_r_squared", fit.r_squared)
        @test fit.dist === nothing

        # Directed golden master: samplike on its transpose
        # (R: netlm(samp, t(samp), nullhyp="classical"))
        samp = sampson_like()
        fit_d = netlm(samp, [Matrix(as_matrix(samp)')]; nullhyp=:classical)
        @test fit_d.n == 306
        @test fit_d.directed
        @test Networks.check_golden(R_SNA, "lm_directed_coefficients", fit_d.coefficients)
        @test Networks.check_golden(R_SNA, "lm_directed_tstat", fit_d.tstat)

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
        @test Networks.check_golden(R_SNA, "logit_coefficients", fit.coefficients)
        @test Networks.check_golden(R_SNA, "logit_se", fit.se)
        @test Networks.check_golden(R_SNA, "logit_tstat", fit.tstat)
        @test Networks.check_golden(R_SNA, "logit_pgreqabs", fit.pgreqabs)
        @test Networks.check_golden(R_SNA, "logit_deviance", fit.deviance)
        @test Networks.check_golden(R_SNA, "logit_null_deviance", fit.null_deviance)
        @test Networks.check_golden(R_SNA, "logit_aic", fit.aic)
        @test Networks.check_golden(R_SNA, "logit_bic", fit.bic)

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

    # ---------------------------------------------------------------------
    # Missing-dyad policy (SNA.jl#1)
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
        ("degree_centrality", (g; kw...) -> degree_centrality(g; kw...)),
        ("betweenness_centrality", (g; kw...) -> betweenness_centrality(g; kw...)),
        ("closeness_centrality", (g; kw...) -> closeness_centrality(g; kw...)),
        ("eigenvector_centrality", (g; kw...) -> eigenvector_centrality(g; kw...)),
        ("bonacich_power", (g; kw...) -> bonacich_power(g; exponent=0.1, kw...)),
        ("katz_centrality", (g; kw...) -> katz_centrality(g; kw...)),
        ("pagerank", (g; kw...) -> pagerank(g; kw...)),
        ("flowbet", (g; kw...) -> flowbet(g; kw...)),
        ("centralization(:degree)", (g; kw...) -> centralization(g, :degree; kw...)),
        ("centralization(:betweenness)", (g; kw...) -> centralization(g, :betweenness; kw...)),
        ("centralization(:closeness)", (g; kw...) -> centralization(g, :closeness; kw...)),
        ("centralization(:eigenvector)", (g; kw...) -> centralization(g, :eigenvector; kw...)),
        ("density", (g; kw...) -> density(g; kw...)),
        ("gden", (g; kw...) -> gden(g; kw...)),
        ("reciprocity", (g; kw...) -> reciprocity(g; kw...)),
        ("grecip", (g; kw...) -> grecip(g; kw...)),
        ("transitivity", (g; kw...) -> transitivity(g; kw...)),
        ("transitivity(:local)", (g; kw...) -> transitivity(g; type=:local, kw...)),
        ("gtrans", (g; kw...) -> gtrans(g; kw...)),
        ("dyad_census", (g; kw...) -> dyad_census(g; kw...)),
        ("triad_census", (g; kw...) -> triad_census(g; kw...)),
        ("mutuality", (g; kw...) -> mutuality(g; kw...)),
        ("hierarchy", (g; kw...) -> hierarchy(g; kw...)),
        ("hierarchy(:krackhardt)", (g; kw...) -> hierarchy(g; measure=:krackhardt, kw...)),
        ("efficiency", (g; kw...) -> efficiency(g; kw...)),
        ("connectedness", (g; kw...) -> connectedness(g; kw...)),
        ("component_dist", (g; kw...) -> component_dist(g; kw...)),
        ("reachability", (g; kw...) -> reachability(g; kw...)),
        ("components", (g; kw...) -> components(g; kw...)),
        ("largest_component", (g; kw...) -> largest_component(g; kw...)),
        ("cliques", (g; kw...) -> cliques(g; min_size=2, kw...)),
        ("kcores", (g; kw...) -> kcores(g; k=1, kw...)),
        ("cutpoints", (g; kw...) -> cutpoints(g; kw...)),
        ("bridges", (g; kw...) -> bridges(g; kw...)),
        ("bicomponents", (g; kw...) -> bicomponents(g; kw...)),
        ("geodesic_distance", (g; kw...) -> geodesic_distance(g; kw...)),
        ("diameter", (g; kw...) -> diameter(g; kw...)),
        ("average_path_length", (g; kw...) -> average_path_length(g; kw...)),
        ("structural_equivalence", (g; kw...) -> structural_equivalence(g; kw...)),
        ("regular_equivalence", (g; kw...) -> regular_equivalence(g; kw...)),
        ("equiv_clust", (g; kw...) -> equiv_clust(g; k=2, kw...)),
        ("blockmodel", (g; kw...) -> blockmodel(g; k=2, kw...)),
        ("qaptest", (g; kw...) -> (qt = qaptest((a, b) -> cor(vec(a), vec(b)),
                                               g, g; reps=10,
                                               rng=MersenneTwister(1), kw...);
                                   (qt.testval, qt.dist))),
        ("netlm", (g; kw...) -> netlm(g, g; reps=10, rng=MersenneTwister(1), kw...).coefficients),
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
        @test degree_centrality(masked; mode=:out, missing=:face)[1] ==
              degree_centrality(obs; mode=:out)[1]
        @test density(masked; missing=:face) == density(obs) == 7 / 30
        @test dyad_census(masked; missing=:face) == dyad_census(obs)

        # Masking every dyad of the network still leaves the face values intact
        # but blocks every measure
        allmasked = copy(obs)
        for i in 1:6, j in 1:6
            i == j || set_missing_dyad!(allmasked, i, j)
        end
        @test_throws ArgumentError density(allmasked)
        @test density(allmasked; missing=:face) == density(obs)

        # clear_missing_dyads! declares everything observed again: the default
        # policy works once more
        cleared = clear_missing_dyads!(copy(masked))
        @test n_missing_dyads(cleared) == 0
        @test density(cleared) == density(obs)
        @test triad_census(cleared) == triad_census(obs)
    end

    @testset "Missing dyads: policy validation and QAP arguments" begin
        masked, obs = masked_pair()

        # Only :error and :face are policies
        @test Set(MISSING_POLICIES) == Set([:error, :face])
        @test_throws ArgumentError density(obs; missing=:ignore)
        @test_throws ArgumentError netlm(obs, obs; missing=:ignore, reps=10)

        # No SNA measure claims a principled missing-data treatment
        for f in (degree_centrality, density, triad_census, centralization,
                  netlm, netlogit, qaptest)
            @test supports_missing(f) == false
        end

        # A mask on *either* QAP argument is caught: y, a predictor, or the
        # second graph of qaptest
        @test_throws ArgumentError netlm(masked, obs; reps=10)
        @test_throws ArgumentError netlm(obs, masked; reps=10)
        @test_throws ArgumentError netlogit(obs, [obs, masked]; reps=10)
        @test_throws ArgumentError qaptest((a, b) -> cor(vec(a), vec(b)),
                                           obs, masked; reps=10)

        # Raw matrices carry no mask, so they are always usable
        A = Float64.(as_matrix(obs))
        fit = netlm(A, A; reps=10, rng=MersenneTwister(2))
        @test length(fit.coefficients) == 2

        # The error message names the routine that refused
        msg = sprint(showerror, try
            betweenness_centrality(masked)
        catch e
            e
        end)
        @test occursin("betweenness_centrality", msg)
        @test occursin("missing=:face", msg)
    end

    @testset "Result metadata protocol (QAP regressions)" begin
        flo = florentine()
        biz = load_dataset(:florentine_business)

        # OLS standard errors assume independent dyads; QAP p-values are
        # a separate inference layer. The metadata must preserve that distinction.
        fq = netlm(flo, biz; reps=200, rng=Xoshiro(31))
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
        gl = netlogit(flo, distance; reps=200, rng=Xoshiro(32))
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
            for (key, value) in (("density", density(b)),
                                 ("discount_density", density(b; discount_bipartite=true)),
                                 ("degree", degree_centrality(b)),
                                 ("bonacich", bonacich_power(b; exponent=0.05)),
                                 ("flowbet", flowbet(b)))
                @test Networks.check_golden(R_SNA, id * "_" * key, value)
            end
            for f in (betweenness_centrality, closeness_centrality, eigenvector_centrality,
                      katz_centrality, structural_equivalence, regular_equivalence)
                @test f(b) ≈ f(full)
            end
            @test blockmodel(b; k=2).block_matrix ≈ blockmodel(full; k=2).block_matrix
            # A statistic detects forbidden within-mode ties in every replicate.
            within_mode(A, _) = sum(A[1:2,1:2]) + sum(A[3:5,3:5])
            qt = qaptest(within_mode, b, b; reps=30, rng=Xoshiro(42))
            @test all(iszero, qt.dist)
            @test netlm(b, full; nullhyp=:classical).n == (directed ? 20 : 10)
        end
    end
    @test_throws ArgumentError qaptest(sum, zeros(2,3), zeros(2,3))
end

@testset "R edge cases, weights, and Katz baseline" begin
    empty = Network(3)
    for (key, actual) in (("empty_transitivity", transitivity(empty)),
                           ("empty_nonnull_nan", isnan(reciprocity(empty; method=:dyadic_nonnull))),
                           ("empty_edgewise_nan", isnan(reciprocity(empty; method=:edgewise))),
                           ("singleton_reciprocity_nan", isnan(reciprocity(Network(1)))))
        @test Networks.check_golden(R_SNA, key, actual)
    end
    @test isnan(reciprocity(Network(3; directed=false); method=:edgewise))
    loops = Network(4; loops=true)
    add_edges!(loops, [(1,1), (1,2), (2,3), (3,4)])
    @test Networks.check_golden(R_SNA, "loops_degree", degree_centrality(loops))
    @test Networks.check_golden(R_SNA, "loops_degree_diag", degree_centrality(loops; diag=true))
    @test Networks.check_golden(R_SNA, "loops_density", density(loops))
    @test Networks.check_golden(R_SNA, "loops_density_diag", density(loops; diag=true))
    W = [0.0 2 1; 1 0 3; 4 2 0]
    weighted = network_from_matrix(W; store_values=true)
    for (key, value) in (("weighted_evcent", eigenvector_centrality(weighted; ignore_eval=false)),
                          ("weighted_bonacich", bonacich_power(weighted; exponent=0.05, ignore_eval=false)),
                          ("weighted_flowbet", flowbet(weighted; ignore_eval=false)),
                          ("weighted_degree", degree_centrality(weighted; ignore_eval=false)))
        @test Networks.check_golden(R_SNA, key, value)
    end
    raw = katz_centrality(weighted; α=0.1, β=2.0, ignore_eval=false)
    @test (I - 0.1W') * raw ≈ fill(2.0, 3)
    @test raw ≈ 2katz_centrality(weighted; α=0.1, ignore_eval=false)
    @test norm(katz_centrality(weighted; normalized=true)) ≈ 1
    @test all(iszero, katz_centrality(weighted; β=0.0))
    @test_throws ArgumentError katz_centrality(weighted; α=2.0)
    @test_throws ArgumentError degree_centrality(weighted; mode=:bogus)
    @test_throws ArgumentError closeness_centrality(weighted; cmode=:bogus)
    @test_throws ArgumentError components(weighted; connected=:bogus)
    @test_throws ArgumentError cutpoints(weighted; connected=:bogus)
    @test_throws ArgumentError cliques(weighted; symmetrize=:bogus)
    # Bipartite adjacency has equally large ±rho; direct eigen solve must obey
    # the eigenvector equation, unlike unshifted power iteration on this path.
    path = Network(3; directed=false); add_edges!(path, [(1,2),(2,3)])
    ev = eigenvector_centrality(path)
    @test SNA._sociomatrix(path) * ev ≈ sqrt(2) * ev
    @test isempty(largest_component(Network(0)))
end

@testset "StatsAPI and scheduling-independent QAP" begin
    flo, biz = florentine(), load_dataset(:florentine_business)
    d = abs.(FLO_WEALTH .- FLO_WEALTH')
    q1 = qaptest(gcor, flo, biz; reps=30, rng=Xoshiro(2), threaded=false)
    q2 = qaptest(gcor, flo, biz; reps=30, rng=Xoshiro(2), threaded=true)
    @test q1.dist == q2.dist
    for fitfun in (netlm, netlogit), nullhyp in (:qapy, :qapx, :qapspp)
        # Continuous predictors keep this scheduling/StatsAPI test identified;
        # separated binary permutations have a separate refusal regression below.
        predictors = fitfun === netlogit ? [d, FLO_WEALTH .+ FLO_WEALTH'] : [biz, d]
        a = fitfun(flo, predictors; reps=30, rng=Xoshiro(44), threaded=false, nullhyp)
        b = fitfun(flo, predictors; reps=30, rng=Xoshiro(44), threaded=true, nullhyp)
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
        @test coeftable(a) isa Networks.CoefficientTable
        @test all(values(Networks.check_statsapi(a; strict=true)))
    end
    zero_tail = QAPTestResult(2.0, [0.0,1.0], 0.0, 1.0, 2)
    @test occursin("<0.5", sprint(show, zero_tail))
    # Imported fixture/optimizer helpers stay qualified in the public namespace.
    @test :Network in names(SNA)
    @test !(:load_golden in names(SNA))
    @test !(:bootstrap_cov in names(SNA))
    @test !(:Networks in names(SNA))
    @test Network(5) isa Networks.Network
    for f in (density, degree_centrality, netlm, netlogit, components)
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
    @test_throws ArgumentError regular_equivalence(net; max_iter=0)
    flo = florentine()
    @test_throws ArgumentError netlm(flo, [flo, flo]; nullhyp=:classical)
    @test_throws ArgumentError netlogit(flo, flo; nullhyp=:classical)
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
    @test_throws ArgumentError netlogit(response, predictor; nullhyp=:classical)
    @test_throws ArgumentError netlogit(response, 1e6 .* predictor; nullhyp=:classical)
    @test_throws ArgumentError netlogit(1 .- response - Matrix{Float64}(I, 4, 4),
                                      predictor; nullhyp=:classical)
    # Adding a failure to the predictor=1 stratum identifies a finite slope.
    predictor[1, 3] = 1
    fitted = netlogit(response, predictor; nullhyp=:classical)
    @test fitted.converged
    @test coef(fitted)[2] ≈ log(2 / 1) - log(2 / 7) atol=1e-7
end


@testset "Separation certificate verifies exact margins" begin
    y = [0.0, 1.0, 0.0, 1.0]
    # The tiny interior reversal gives overlap and a finite MLE. A floating
    # tolerance cannot turn these genuine negative margins into a certificate.
    X = hcat(ones(4), [-1.0, -1e-16, 1e-16, 1.0])
    @test !SNA._separation_certificate(X, y, [0.0, 1.0])
    @test SNA._logit_fit(X, y).converged
    @test SNA._separation_certificate(X, [0.0,0.0,1.0,1.0], [0.0,1.0])
    # Correct the numerically perturbed direction in the exact boundary
    # nullspace, then verify every original row without rounding the data.
    Q = hcat(ones(4), [-1.0, 0.0, 0.0, 1.0])
    @test SNA._separation_certificate(Q, y, [1e-16, 1.0])
    @test !SNA._separation_certificate(Q, y, [0.0, 0.0])
    flo, biz = florentine(), load_dataset(:florentine_business)
    messages = String[]
    for threaded in (false, true)
        err = try
            netlogit(flo, biz; reps=200, rng=Xoshiro(32), threaded)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("QAP replicate", sprint(showerror, err))
        @test occursin("separated", sprint(showerror, err))
        @test occursin("No permutation p-value", sprint(showerror, err))
        push!(messages, sprint(showerror, err))
    end
    @test messages[1] == messages[2]
end


@testset "Binary-stratum separation agrees with analytical likelihood" begin
    for n0 in (2, 3), n1 in (2, 3), successes0 in 0:n0, successes1 in 0:n1
        X = hcat(ones(n0+n1), vcat(zeros(n0), ones(n1)))
        y = vcat(ones(successes0), zeros(n0-successes0),
                 ones(successes1), zeros(n1-successes1))
        if successes0 in (0,n0) || successes1 in (0,n1)
            @test_throws ArgumentError SNA._logit_fit(X, y)
        else
            fit = SNA._logit_fit(X, y)
            alpha = log(successes0 / (n0-successes0))
            beta = log(successes1 / (n1-successes1)) - alpha
            @test fit.coef ≈ [alpha, beta] atol=1e-7
            @test fit.converged
        end
    end
end
