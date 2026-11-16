# Trees from distance matrices: UPGMA, WPGMA, neighbour joining, BIONJ.

# Saturated pairs (Inf/NaN) can't be joined; give them a distance a little
# beyond the largest finite one, with a warning.
function _finite_distances(d::DistanceMatrix)
    D = copy(d.d)
    bad = .!isfinite.(D)
    for i in axes(D, 1)
        bad[i, i] = false
        D[i, i] = 0.0
    end
    if any(bad)
        fin = filter(isfinite, D)
        big = 2 * maximum(fin; init=1.0)
        D[bad] .= big
        @warn "$(count(bad) ÷ 2) pairs are saturated or share no sites; their distance was set to $(round(big; sigdigits=3)). Consider model = :p or a gamma model."
    end
    return D
end

"""
    upgma(d; root=nothing, ladderize=false) -> PhyloTree
    upgma(aln; model=:JC69, gamma=nothing, kwargs...) -> PhyloTree

UPGMA (average linkage) tree from a [`DistanceMatrix`](@ref), or from an
alignment via [`distances`](@ref). The tree is rooted and ultrametric.
"""
upgma(d::DistanceMatrix; kwargs...) = _cluster(d, :upgma; kwargs...)

"""
    wpgma(d; root=nothing, ladderize=false) -> PhyloTree

WPGMA (McQuitty) tree: as [`upgma`](@ref), but a merged cluster's distance
is the simple mean of its two parts' distances.
"""
wpgma(d::DistanceMatrix; kwargs...) = _cluster(d, :wpgma; kwargs...)

function _cluster(d::DistanceMatrix, method::Symbol; root=nothing, ladderize::Bool=false)
    n = ntaxa(d)
    n == 0 && throw(ArgumentError("no taxa"))
    D = _finite_distances(d)
    total = 2n - 1
    children = [Int[] for _ in 1:total]
    height = zeros(total)
    len = fill(NaN, total)
    size_ = ones(Int, total)
    alive = trues(n)
    node = collect(1:n)                 # tree node of each live row
    # nearest later neighbour of each row (first on ties), so a merge only
    # rescans the rows that pointed at the merged pair
    nn = zeros(Int, n)
    nnd = fill(Inf, n)
    foreach(i -> _rescan!(nn, nnd, D, alive, i), 1:n)
    nxt = n
    for _ in 1:n-1
        # closest pair (first in row order on ties, as R's hclust)
        i = 0
        for k in 1:n
            alive[k] && nn[k] != 0 && (i == 0 || nnd[k] < nnd[i]) && (i = k)
        end
        j = nn[i]
        best = D[i, j]
        u = (nxt += 1)
        h = best / 2
        ni, nj = node[i], node[j]
        children[u] = [ni, nj]
        height[u] = h
        len[ni] = h - height[ni]
        len[nj] = h - height[nj]
        size_[u] = size_[ni] + size_[nj]
        alive[j] = false
        for k in 1:n
            (alive[k] && k != i) || continue
            v = method === :upgma ? (size_[ni] * D[i, k] + size_[nj] * D[j, k]) / size_[u] :
                (D[i, k] + D[j, k]) / 2
            D[i, k] = D[k, i] = v
        end
        node[i] = u
        for k in 1:n
            alive[k] || continue
            if k == i || nn[k] == i || nn[k] == j
                _rescan!(nn, nnd, D, alive, k)
            elseif k < i && (D[k, i] < nnd[k] || (D[k, i] == nnd[k] && i < nn[k]))
                nn[k] = i; nnd[k] = D[k, i]
            end
        end
    end
    labels = vcat(d.names, fill("", total - n))
    t = _build(children, len, labels, total)
    return _finish(t; root, ladderize)
end

# the closest live row after `i` (the first on ties), or 0
function _rescan!(nn, nnd, D, alive, i)
    nn[i] = 0; nnd[i] = Inf
    for k in i+1:length(alive)
        if alive[k] && D[i, k] < nnd[i]
            nn[i] = k; nnd[i] = D[i, k]
        end
    end
end

"""
    nj(d; root=nothing, nonnegative=true, ladderize=false) -> PhyloTree
    nj(aln; model=:JC69, gamma=nothing, kwargs...) -> PhyloTree

Neighbour-joining tree (Saitou & Nei 1987, in Studier & Keppler's form)
from a [`DistanceMatrix`](@ref) or an alignment. The tree is unrooted
(drawn from a node with three children) unless `root` is `:midpoint` or an
outgroup. With `nonnegative = true` negative branch lengths are set to 0.
"""
nj(d::DistanceMatrix; kwargs...) = _nj(d, false; kwargs...)

"""
    bionj(d; root=nothing, nonnegative=true, ladderize=false) -> PhyloTree
    bionj(aln; model=:JC69, gamma=nothing, kwargs...) -> PhyloTree

BIONJ tree (Gascuel 1997): neighbour joining that weights each merge by
the estimated variance of the distances, which is more accurate when
distances are large. Same keywords as [`nj`](@ref).
"""
bionj(d::DistanceMatrix; kwargs...) = _nj(d, true; kwargs...)

function _nj(d::DistanceMatrix, bio::Bool; root=nothing, nonnegative::Bool=true, ladderize::Bool=false)
    n = ntaxa(d)
    n == 0 && throw(ArgumentError("no taxa"))
    if n <= 2
        children = n == 1 ? [Int[]] : [Int[], Int[], [1, 2]]
        len = n == 1 ? [NaN] : (l = d.d[1, 2] / 2; [l, l, NaN])
        return _finish(_build(children, len, vcat(d.names, fill("", length(children) - n)), length(children)); root, ladderize)
    end
    D = _finite_distances(d)
    V = bio ? copy(D) : D
    total = 2n - 2
    children = [Int[] for _ in 1:total]
    len = fill(NaN, total)
    active = collect(1:n)
    node = collect(1:n)
    S = vec(sum(D; dims=2))
    nxt = n
    while length(active) > 3
        r = length(active)
        best = Inf; bi = bj = 0
        @inbounds for a in 1:r
            i = active[a]
            for b in a+1:r
                j = active[b]
                q = (r - 2) * D[i, j] - S[i] - S[j]
                if q < best
                    best = q; bi = a; bj = b
                end
            end
        end
        i, j = active[bi], active[bj]
        dij = D[i, j]
        li = (dij + (S[i] - S[j]) / (r - 2)) / 2
        lj = dij - li
        λ = 0.5
        if bio && V[i, j] > 0
            s = 0.0
            for k in active
                (k == i || k == j) && continue
                s += V[j, k] - V[i, k]
            end
            λ = clamp(0.5 + s / (2 * (r - 2) * V[i, j]), 0.0, 1.0)
        end
        u = (nxt += 1)
        children[u] = [node[i], node[j]]
        len[node[i]] = li
        len[node[j]] = lj
        # row i becomes the new node
        for k in active
            (k == i || k == j) && continue
            dik, djk = D[i, k], D[j, k]
            if bio
                nd = λ * (dik - li) + (1 - λ) * (djk - lj)
                V[i, k] = V[k, i] = λ * V[i, k] + (1 - λ) * V[j, k] - λ * (1 - λ) * V[i, j]
            else
                nd = (dik + djk - dij) / 2
            end
            S[k] += nd - dik - djk
            D[i, k] = D[k, i] = nd
        end
        node[i] = u
        deleteat!(active, bj)
        S[i] = _rowsum(D, i, active)
    end
    a, b, c = active
    u = (nxt += 1)
    children[u] = [node[a], node[b], node[c]]
    len[node[a]] = (D[a, b] + D[a, c] - D[b, c]) / 2
    len[node[b]] = (D[a, b] + D[b, c] - D[a, c]) / 2
    len[node[c]] = (D[a, c] + D[b, c] - D[a, b]) / 2
    nonnegative && (len .= ifelse.(len .< 0, 0.0, len))
    labels = vcat(d.names, fill("", total - n))
    t = _build(children, len, labels, u)
    return _finish(t; root, ladderize)
end

# (a function rather than a generator inside `_nj`: a closure capturing the
# reassigned `i` would box it and make the whole search type-unstable)
function _rowsum(D, i, active)
    s = 0.0
    for k in active
        k == i || (s += D[i, k])
    end
    return s
end

# the alignment forms: distances first
for f in (:upgma, :wpgma, :nj, :bionj)
    @eval function $f(a::Alignment; model::Symbol=:JC69, gamma=nothing,
                      progress::Bool=_progress_default(), kwargs...)
        return $f(distances(a; model, gamma, progress); kwargs...)
    end
end
