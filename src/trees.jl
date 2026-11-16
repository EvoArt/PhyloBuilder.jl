# Building PhyloIO trees, and rooting them.

# A PhyloTree from child lists: node `root` at the top, `len[i]` the branch
# above node `i`, `label[i]` its label (tips are usually 1:n).
function _build(children::AbstractVector, len::AbstractVector, label::AbstractVector, root::Integer;
                data=nothing)
    t = PhyloTree()
    stack = [(Int(root), 0)]
    while !isempty(stack)
        i, p = pop!(stack)
        j = addnode!(t; parent=p, length=p == 0 ? NaN : len[i], label=label[i])
        data === nothing || isempty(data[i]) || merge!(t.data[j], data[i])
        ch = children[i]
        for k in length(ch):-1:1
            push!(stack, (ch[k], j))
        end
    end
    return t
end

# ------------------------------------------------------------------ unrooted views

# The tree as an undirected graph: neighbours (children first, in drawing
# order, then the parent), branch lengths and edge annotations (stored at the
# child end in PhyloTree, here per edge). A root with two children is
# dissolved into one edge, so the graph is the unrooted tree.
struct _Graph
    nbr::Vector{Vector{Int}}
    len::Dict{Tuple{Int,Int},Float64}
    data::Dict{Tuple{Int,Int},Dict{Symbol,Any}}
    label::Vector{String}
    nodedata::Vector{Dict{Symbol,Any}}
    gone::Int                      # a dissolved root, or 0
end

_ekey(a, b) = a < b ? (a, b) : (b, a)
_elen(g::_Graph, a, b) = g.len[_ekey(a, b)]

function _graph(t::PhyloTree)
    n = nnodes(t)
    nbr = [Int[] for _ in 1:n]
    len = Dict{Tuple{Int,Int},Float64}()
    data = Dict{Tuple{Int,Int},Dict{Symbol,Any}}()
    for i in 1:n
        append!(nbr[i], childnodes(t, i))
        p = parentnode(t, i)
        p == 0 && continue
        push!(nbr[i], p)
        l = branchlength(t, i)
        len[_ekey(i, p)] = isnan(l) ? 0.0 : l
        data[_ekey(i, p)] = copy(t.data[i])
    end
    r = rootnode(t)
    gone = 0
    if length(childnodes(t, r)) == 2
        a, b = childnodes(t, r)
        nbr[a][findfirst(==(r), nbr[a])] = b
        nbr[b][findfirst(==(r), nbr[b])] = a
        len[_ekey(a, b)] = len[_ekey(a, r)] + len[_ekey(b, r)]
        da, db = data[_ekey(a, r)], data[_ekey(b, r)]
        data[_ekey(a, b)] = isleaf(t, a) ? db : da
        delete!(len, _ekey(a, r)); delete!(len, _ekey(b, r))
        empty!(nbr[r])
        gone = r
    end
    nodedata = [Dict{Symbol,Any}() for _ in 1:n]   # node-only annotations: none kept
    return _Graph(nbr, len, data, [nodelabel(t, i) for i in 1:n], nodedata, gone)
end

# A PhyloTree hanging from node `r` of the graph (r ≠ the dissolved root).
function _tree_from(g::_Graph, r::Int)
    t = PhyloTree()
    stack = [(r, 0, 0)]                        # node, graph parent, tree parent
    while !isempty(stack)
        i, from, p = pop!(stack)
        j = addnode!(t; parent=p, length=from == 0 ? NaN : _elen(g, i, from), label=g.label[i])
        from == 0 || merge!(t.data[j], g.data[_ekey(i, from)])
        ch = filter(!=(from), g.nbr[i])
        for k in length(ch):-1:1
            push!(stack, (ch[k], i, j))
        end
    end
    return t
end

# A PhyloTree rooted on edge u–v, at distance `x` from `u`.
function _tree_on_edge(g::_Graph, u::Int, v::Int, x::Real)
    L = _elen(g, u, v)
    d = g.data[_ekey(u, v)]
    t = PhyloTree()
    root = addnode!(t)
    for (a, b, l) in ((u, v, x), (v, u, L - x))
        stack = [(a, b, root, Float64(l))]
        while !isempty(stack)
            i, from, p, li = pop!(stack)
            j = addnode!(t; parent=p, length=li, label=g.label[i])
            merge!(t.data[j], p == root ? d : g.data[_ekey(i, from)])
            ch = filter(!=(from), g.nbr[i])
            for k in length(ch):-1:1
                push!(stack, (ch[k], i, j, _elen(g, ch[k], i)))
            end
        end
    end
    return t
end

# distances from `s` to every node, and the predecessor on the path
function _paths(g::_Graph, s::Int)
    n = length(g.nbr)
    dist = fill(Inf, n)
    pred = zeros(Int, n)
    dist[s] = 0
    stack = [s]
    while !isempty(stack)
        i = pop!(stack)
        for j in g.nbr[i]
            isfinite(dist[j]) && continue
            dist[j] = dist[i] + _elen(g, i, j)
            pred[j] = i
            push!(stack, j)
        end
    end
    return dist, pred
end

"""
    midpoint(tree) -> PhyloTree

Root `tree` at the midpoint of its longest tip-to-tip path. Support values
and other branch annotations stay with their branches.
"""
function midpoint(t::PhyloTree)
    ntips(t) < 2 && return copy(t)
    g = _graph(t)
    tp = tips(t)
    d0, _ = _paths(g, tp[1])
    a = tp[argmax([d0[i] for i in tp])]
    da, pred = _paths(g, a)
    b = tp[argmax([da[i] for i in tp])]
    half = da[b] / 2
    # walk from b towards a until the path passes the midpoint
    u = b
    while da[pred[u]] > half
        u = pred[u]
    end
    v = pred[u]
    # the midpoint lies on edge u–v at distance da[u] - half from u
    return _tree_on_edge(g, u, v, da[u] - half)
end

"""
    unroot(tree) -> PhyloTree

`tree` with its root dissolved: a root with two children is replaced by a
single branch, so the tree is drawn from a node with three or more
children, as unrooted trees are written in Newick.
"""
function unroot(t::PhyloTree)
    length(childnodes(t, rootnode(t))) == 2 || return copy(t)
    g = _graph(t)
    # hang it from the first internal node, in pre-order
    r = first(i for i in preorder(t) if i != g.gone && !isleaf(t, i))
    return _tree_from(g, r)
end

"""
    root(tree, outgroup) -> PhyloTree
    root(tree, :midpoint) -> PhyloTree

Root `tree` on the branch above `outgroup` (a tip label, or several whose
common ancestor is used, halving that branch), or at the midpoint.
"""
root(t::PhyloTree, how::Symbol) = how === :midpoint ? midpoint(t) :
    throw(ArgumentError("root: :midpoint or an outgroup"))
function root(t::PhyloTree, outgroup)
    o = nodeindex(t, outgroup)
    if !(outgroup isa Integer) && !(outgroup isa AbstractString) && length(outgroup) > 1
        # the MRCA in a tree rooted inside the outgroup is not the outgroup's
        # clade: root it on a tip outside first
        names = Set(String.(outgroup))
        other = first(i for i in tips(t) if !(nodelabel(t, i) in names))
        t = reroot(t, other)
        o = mrca(t, collect(names))
    end
    return reroot(t, o)
end

# the `root`/`ladderize` keywords of the tree builders
function _finish(t::PhyloTree; root=nothing, ladderize::Bool=false)
    if root === :midpoint
        t = midpoint(t)
    elseif root !== nothing && root !== :none
        t = PhyloBuilder.root(t, root)
    end
    ladderize && ladderize!(t)
    return t
end
