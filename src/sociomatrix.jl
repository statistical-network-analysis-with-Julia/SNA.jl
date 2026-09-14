# Algorithms operate on all actors in both modes, never on an incidence matrix.
_graph(net::Network) = net.graph
_graph(net::BipartiteNetwork) = net.network.graph

function _sociomatrix(net::AbstractNetwork; ignore_eval::Bool=true,
                     attr::Symbol=:weight, diag::Bool=false)
    A = Matrix{Float64}(as_matrix(net; expand_bipartite=true,
                                 attr=ignore_eval ? nothing : attr))
    all(isfinite, A) || throw(ArgumentError("edge attribute $attr must contain finite values"))
    if !diag
        for i in axes(A, 1)
            A[i, i] = 0.0
        end
    end
    return A
end

function _sociomatrix(x::AbstractMatrix)
    size(x, 1) == size(x, 2) ||
        throw(ArgumentError("a sociomatrix must be square; expand two-mode data to all actors"))
    A = Matrix{Float64}(x)
    all(isfinite, A) || throw(ArgumentError("a sociomatrix must contain finite values"))
    return A
end

function _symmetrized_graph(net, rule::Symbol)
    rule in (:strong, :weak) ||
        throw(ArgumentError("symmetrize must be :strong (mutual arcs) or :weak (either arc)"))
    g = Graphs.SimpleGraph(nv(net))
    for e in edges(net)
        i, j = src(e), dst(e)
        i == j && continue
        if !is_directed(net) || rule == :weak || has_edge(net, j, i)
            Graphs.add_edge!(g, i, j)
        end
    end
    return g
end
