# Pairwise evolutionary distances between aligned sequences.

"""
    DistanceMatrix(names, d::AbstractMatrix; model=:unknown)

A symmetric matrix of pairwise distances with the names of the sequences,
as returned by [`distances`](@ref). It is an `AbstractMatrix{Float64}`, so
it can be used anywhere a matrix can; `d["a", "b"]` looks up a pair by name.
Pass it to [`nj`](@ref), [`bionj`](@ref), [`upgma`](@ref) or
[`wpgma`](@ref) to build a tree.
"""
struct DistanceMatrix <: AbstractMatrix{Float64}
    names::Vector{String}
    d::Matrix{Float64}
    model::Symbol
    gamma::Union{Nothing,Float64}
    function DistanceMatrix(names::AbstractVector, d::AbstractMatrix; model::Symbol=:unknown, gamma=nothing)
        n = length(names)
        size(d) == (n, n) || throw(ArgumentError("a $(size(d)) matrix for $n names"))
        allunique(names) || throw(ArgumentError("names must be unique"))
        new(String.(string.(names)), Matrix{Float64}(d), model, gamma === nothing ? nothing : Float64(gamma))
    end
end

Base.size(d::DistanceMatrix) = size(d.d)
Base.getindex(d::DistanceMatrix, i::Int, j::Int) = d.d[i, j]
Base.getindex(d::DistanceMatrix, a::AbstractString, b::AbstractString) = d.d[_dindex(d, a), _dindex(d, b)]
Base.IndexStyle(::Type{DistanceMatrix}) = IndexCartesian()
Base.copy(d::DistanceMatrix) = DistanceMatrix(d.names, copy(d.d); model=d.model, gamma=d.gamma)
taxa(d::DistanceMatrix) = d.names
ntaxa(d::DistanceMatrix) = length(d.names)

function _dindex(d::DistanceMatrix, a::AbstractString)
    i = findfirst(==(a), d.names)
    i === nothing && throw(KeyError(a))
    return i
end

const _DIST_MODELS = (:raw, :p, :JC69, :K80, :F81, :TN93, :logdet, :paralin)

"""
    distances(aln; model=:JC69, gamma=nothing, threaded=true, progress) -> DistanceMatrix

Pairwise distances between the sequences of `aln`. Sites where either
sequence has a gap or an ambiguous base are left out of that pair
(pairwise deletion).

`model` is one of

- `:raw` – number of differences
- `:p` – proportion of sites that differ
- `:JC69` – Jukes & Cantor (1969)
- `:K80` – Kimura (1980) two-parameter, transitions vs transversions
- `:F81` – Felsenstein (1981), unequal base frequencies (A, C, G and T
  counted over the whole alignment, ignoring ambiguity codes)
- `:TN93` – Tamura & Nei (1993), two transition rates and unequal frequencies
- `:logdet` – LogDet (Lockhart et al. 1994; Gu & Li 1996), robust to unequal composition
- `:paralin` – paralinear (Lake 1994), LogDet corrected for the composition of each sequence

`gamma = α` adds gamma-distributed rates across sites with shape `α`
(JC69, K80, F81, TN93). Pairs too divergent for the model (saturated)
get `Inf`. The results match R's `ape::dist.dna` (`pairwise.deletion = TRUE`).
"""
function distances(a::Alignment; model::Symbol=:JC69, gamma=nothing, threaded::Bool=true,
                   progress::Bool=_progress_default())
    freqs = model in (:F81, :TN93) ? basefreqs(a; ambiguous=false) : nothing
    return _distances(a.names, _patterns(a), freqs; model, gamma, threaded,
                      prog=_progress(ntaxa(a), "Distances"; unit="rows", enabled=progress && ntaxa(a) >= 200))
end

function _distances(names, pat::_Patterns, freqs; model::Symbol=:JC69, gamma=nothing, threaded::Bool=true,
                    prog=_NOPROGRESS)
    model in _DIST_MODELS || throw(ArgumentError("model must be one of $(join(_DIST_MODELS, ", "))"))
    gamma === nothing || model in (:JC69, :K80, :F81, :TN93) ||
        throw(ArgumentError("gamma is available for JC69, K80, F81 and TN93"))
    n = length(names)
    bb = _basebits(pat)
    D = zeros(n, n)
    comp = model === :paralin ? _composition(bb) : nothing
    rowwork = function (i)
        N = zeros(4, 4)
        for j in i+1:n
            _paircounts!(N, bb, i, j)
            marg = comp === nothing ? nothing : (view(comp, :, i), view(comp, :, j))
            D[i, j] = D[j, i] = _pairdist(N, model, gamma, freqs, marg)
        end
        _next!(prog)
    end
    if threaded && Threads.nthreads() > 1 && n > 20
        Threads.@threads :dynamic for i in 1:n
            rowwork(i)
        end
    else
        foreach(rowwork, 1:n)
    end
    nsat = count(!isfinite, D) ÷ 2
    _finish!(prog; info=nsat > 0 ? _sty(prog.io, "$nsat saturated pairs"; color=_PAL.warn) : nothing)
    return DistanceMatrix(names, D; model, gamma)
end

# The unambiguous bases of every sequence as bit masks: bit s of
# `bits[q, b, i]` says that taxon i has base b at the s-th pattern of word q.
# Gaps and ambiguity codes set no bit, so they drop out of every count
# (pairwise deletion). Patterns are grouped by weight, each group starting on
# a fresh word, so the base-pair counts of a pair are weighted popcounts:
# 64 patterns per instruction instead of one.
struct _BaseBits
    bits::Array{UInt64,3}                          # words × 4 × ntaxa
    groups::Vector{Tuple{UnitRange{Int},Float64}}  # word range and weight
end

function _basebits(pat::_Patterns)
    n = size(pat.data, 1)
    byweight = Dict{Float64,Vector{Int}}()
    for (k, w) in pairs(pat.weight)
        w == 0 || push!(get!(byweight, w, Int[]), k)
    end
    groups = Tuple{UnitRange{Int},Float64}[]
    members = Vector{Int}[]
    nwords = 0
    for w in sort!(collect(keys(byweight)))
        ks = byweight[w]
        push!(groups, (nwords+1:nwords+cld(length(ks), 64), w))
        push!(members, ks)
        nwords += cld(length(ks), 64)
    end
    bits = zeros(UInt64, nwords, 4, n)
    for ((words, _), ks) in zip(groups, members), (s, k) in enumerate(ks)
        q = first(words) + (s - 1) >> 6
        bit = UInt64(1) << ((s - 1) & 63)
        for i in 1:n
            x = pat.data[i, k]
            _isbase(x) && (bits[q, _index(x), i] |= bit)
        end
    end
    return _BaseBits(bits, groups)
end

# N[a, b] = weight of the sites where taxon i has base a and taxon j base b
function _paircounts!(N::Matrix{Float64}, bb::_BaseBits, i::Int, j::Int)
    fill!(N, 0.0)
    B = bb.bits
    for (words, w) in bb.groups
        acc = ntuple(_ -> 0, Val(16))
        @inbounds for q in words
            x = (B[q, 1, i], B[q, 2, i], B[q, 3, i], B[q, 4, i])
            y = (B[q, 1, j], B[q, 2, j], B[q, 3, j], B[q, 4, j])
            # cell k is N[a, b] in column-major order: a = (k-1) % 4 + 1, b = (k-1) ÷ 4 + 1
            acc = acc .+ ntuple(k -> count_ones(x[(k-1)&3+1] & y[(k-1)>>2+1]), Val(16))
        end
        for k in 1:16
            N[k] += w * acc[k]
        end
    end
    return N
end

# base composition of each sequence over its unambiguous sites (4 × ntaxa)
function _composition(bb::_BaseBits)
    n = size(bb.bits, 3)
    c = zeros(4, n)
    for i in 1:n, b in 1:4, (words, w) in bb.groups
        c[b, i] += w * sum(q -> count_ones(bb.bits[q, b, i]), words; init=0)
    end
    return c ./ sum(c; dims=1)
end

# determinant of a 4×4 matrix by cofactors of 2×2 minors (no allocation,
# unlike LinearAlgebra.det's LU)
function _det4(m::AbstractMatrix)
    @inbounds begin
        s0 = m[1, 1] * m[2, 2] - m[2, 1] * m[1, 2]
        s1 = m[1, 1] * m[2, 3] - m[2, 1] * m[1, 3]
        s2 = m[1, 1] * m[2, 4] - m[2, 1] * m[1, 4]
        s3 = m[1, 2] * m[2, 3] - m[2, 2] * m[1, 3]
        s4 = m[1, 2] * m[2, 4] - m[2, 2] * m[1, 4]
        s5 = m[1, 3] * m[2, 4] - m[2, 3] * m[1, 4]
        c5 = m[3, 3] * m[4, 4] - m[4, 3] * m[3, 4]
        c4 = m[3, 2] * m[4, 4] - m[4, 2] * m[3, 4]
        c3 = m[3, 2] * m[4, 3] - m[4, 2] * m[3, 3]
        c2 = m[3, 1] * m[4, 4] - m[4, 1] * m[3, 4]
        c1 = m[3, 1] * m[4, 3] - m[4, 1] * m[3, 3]
        c0 = m[3, 1] * m[4, 2] - m[4, 1] * m[3, 2]
    end
    return s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0
end

_fin(x) = isfinite(x) && x >= 0 ? x : x > -1e-12 ? 0.0 : Inf

function _pairdist(N::AbstractMatrix, model::Symbol, gamma, freqs, marg=nothing)
    L = sum(N)
    L == 0 && return NaN
    same = N[1, 1] + N[2, 2] + N[3, 3] + N[4, 4]
    diffs = L - same
    model === :raw && return diffs
    p = diffs / L
    model === :p && return p
    α = gamma
    # A=1 C=2 G=3 T=4; transitions are A<->G and C<->T
    if model === :JC69
        w = 1 - 4p / 3
        w <= 0 && return Inf
        return _fin(α === nothing ? -0.75 * log(w) : 0.75α * (w^(-1 / α) - 1))
    elseif model === :K80
        P = (N[1, 3] + N[3, 1] + N[2, 4] + N[4, 2]) / L
        Q = p - P
        w1 = 1 - 2P - Q; w2 = 1 - 2Q
        (w1 <= 0 || w2 <= 0) && return Inf
        return _fin(α === nothing ? -0.5log(w1) - 0.25log(w2) :
                    0.5α * (w1^(-1 / α) - 1) + 0.25α * (w2^(-1 / α) - 1))
    elseif model === :F81
        B = 1 - sum(abs2, freqs)
        w = 1 - p / B
        w <= 0 && return Inf
        return _fin(α === nothing ? -B * log(w) : B * α * (w^(-1 / α) - 1))
    elseif model === :TN93
        πA, πC, πG, πT = freqs
        πR = πA + πG; πY = πC + πT
        P1 = (N[1, 3] + N[3, 1]) / L
        P2 = (N[2, 4] + N[4, 2]) / L
        Q = p - P1 - P2
        k1 = 2πA * πG / πR
        k2 = 2πC * πT / πY
        k3 = 2 * (πR * πY - πA * πG * πY / πR - πC * πT * πR / πY)
        w1 = 1 - P1 / k1 - Q / (2πR)
        w2 = 1 - P2 / k2 - Q / (2πY)
        w3 = 1 - Q / (2πR * πY)
        (w1 <= 0 || w2 <= 0 || w3 <= 0) && return Inf
        if α === nothing
            return _fin(-k1 * log(w1) - k2 * log(w2) - k3 * log(w3))
        end
        k4 = 2 * (πA * πG + πC * πT + πR * πY)
        b = -1 / α
        return _fin(α * (k1 * w1^b + k2 * w2^b + k3 * w3^b - k4))
    end
    # LogDet family; F is the divergence matrix of base pairs. As in ape,
    # these are not clamped at zero.
    dN = _det4(N)
    dN <= 0 && return Inf
    logdetF = log(dN) - 4log(L)         # F = N / L
    model === :logdet && return -0.25 * logdetF - log(4)
    # :paralin, as ape: -(log det F - (log det Πx + log det Πy) / 2) / 4, where
    # Πx and Πy are the compositions of the two whole sequences
    fx, fy = marg
    (any(<=(0), fx) || any(<=(0), fy)) && return Inf
    return -0.25 * (logdetF - 0.5 * (sum(log, fx) + sum(log, fy)))
end

"""
    patristic(tree) -> DistanceMatrix

Path-length distances between the tips of a `PhyloTree` (the sum of the
branch lengths between each pair; missing lengths count as 0, or 1 when
the tree has none).
"""
function patristic(t::PhyloTree)
    tp = tips(t)
    n = length(tp)
    haslen = any(i -> !isroot(t, i) && !isnan(branchlength(t, i)), 1:nnodes(t))
    len(i) = haslen ? (isnan(branchlength(t, i)) ? 0.0 : branchlength(t, i)) : 1.0
    depth = rootdistance(t)
    D = zeros(n, n)
    for a in 1:n, b in a+1:n
        m = mrca(t, [tp[a], tp[b]])
        D[a, b] = D[b, a] = depth[tp[a]] + depth[tp[b]] - 2depth[m]
    end
    return DistanceMatrix([nodelabel(t, i) for i in tp], D; model=:patristic)
end
