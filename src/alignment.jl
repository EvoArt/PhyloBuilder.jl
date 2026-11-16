# Aligned nucleotide sequences, stored as 4-bit state sets.
#
# Each site of each sequence is a `UInt8` whose low four bits are the set of
# bases it may be (A = 1, C = 2, G = 4, T/U = 8), so IUPAC ambiguity codes are
# unions (R = A|G = 5, N = 15). A gap is 0x1f: "any base" with a fifth bit set
# so that it can be printed as '-' again. Every computation uses `x & 0x0f`,
# treating gaps as missing data, as PAUP*, PhyML, RAxML and phangorn do.

const _GAP = 0x1f
const _ENC = let t = fill(0xff, 256)
    for (c, v) in ('A' => 0x01, 'C' => 0x02, 'G' => 0x04, 'T' => 0x08, 'U' => 0x08,
                   'R' => 0x05, 'Y' => 0x0a, 'S' => 0x06, 'W' => 0x09, 'K' => 0x0c, 'M' => 0x03,
                   'B' => 0x0e, 'D' => 0x0d, 'H' => 0x0b, 'V' => 0x07, 'N' => 0x0f, 'X' => 0x0f,
                   '?' => 0x0f, 'O' => 0x0f, '-' => _GAP, '.' => _GAP, '~' => _GAP)
        t[UInt8(c)+1] = v
        t[UInt8(lowercase(c))+1] = v
    end
    t
end
const _DEC = let t = fill('N', 32)
    for (v, c) in (0x01 => 'A', 0x02 => 'C', 0x04 => 'G', 0x08 => 'T', 0x05 => 'R', 0x0a => 'Y',
                   0x06 => 'S', 0x09 => 'W', 0x0c => 'K', 0x03 => 'M', 0x0e => 'B', 0x0d => 'D',
                   0x0b => 'H', 0x07 => 'V', 0x0f => 'N', _GAP => '-')
        t[v+1] = c
    end
    t
end

_isbase(x::UInt8) = count_ones(x & 0x0f) == 1
_index(x::UInt8) = trailing_zeros(x) + 1          # A,C,G,T -> 1..4 for a single base

"""
    Alignment(names, sequences)
    Alignment(sequences)                     # names "seq1", "seq2", …
    Alignment(["name" => seq, …])

Aligned nucleotide sequences (all the same length). Sequences can be
strings or anything whose `string` is the sequence (e.g. BioSequences'
`LongDNA`). IUPAC ambiguity codes are kept as sets of bases; gaps (`-`, `.`)
and `N`/`?` are missing data.

`aln[i]` is the `i`th sequence as a `String`, `aln["name"]` looks one up by
name, `aln[taxa, sites]` subsets (by indices, names or masks). Read one with
[`readalignment`](@ref), write one with [`writefasta`](@ref) or
[`writephylip`](@ref).
"""
struct Alignment
    names::Vector{String}
    data::Matrix{UInt8}          # ntaxa × nsites
    function Alignment(names::Vector{String}, data::Matrix{UInt8})
        length(names) == size(data, 1) || throw(ArgumentError("$(length(names)) names for $(size(data, 1)) sequences"))
        allunique(names) || throw(ArgumentError("sequence names must be unique; repeated: " *
                                                join(unique(filter(n -> count(==(n), names) > 1, names)), ", ")))
        new(names, data)
    end
end

function Alignment(names::AbstractVector, seqs::AbstractVector)
    length(names) == length(seqs) || throw(ArgumentError("$(length(names)) names for $(length(seqs)) sequences"))
    strs = [String(string(s)) for s in seqs]
    isempty(strs) && return Alignment(String[], zeros(UInt8, 0, 0))
    lens = map(s -> ncodeunits(filter(!isspace, s)), strs)
    if !allequal(lens)
        lo, hi = extrema(lens)
        throw(ArgumentError("sequences are not aligned (lengths $(lo)–$(hi)); align them first with `align`"))
    end
    data = Matrix{UInt8}(undef, length(strs), lens[1])
    for (i, s) in enumerate(strs)
        j = 0
        for b in codeunits(s)
            isspace(Char(b)) && continue
            v = _ENC[b+1]
            v == 0xff && throw(ArgumentError("sequence $(names[i]): '$(Char(b))' is not a nucleotide code"))
            data[i, j+=1] = v
        end
    end
    return Alignment(String.(string.(names)), data)
end
Alignment(seqs::AbstractVector) = Alignment(["seq$i" for i in eachindex(seqs)], seqs)
Alignment(pairs::AbstractVector{<:Pair}) = Alignment(first.(pairs), last.(pairs))
Alignment(d::AbstractDict) = Alignment(collect(keys(d)), collect(values(d)))

"""    ntaxa(aln)  number of sequences."""
ntaxa(a::Alignment) = size(a.data, 1)
"""    nsites(aln)  alignment length (columns)."""
nsites(a::Alignment) = size(a.data, 2)
"""    taxa(aln)  sequence names."""
taxa(a::Alignment) = a.names

Base.length(a::Alignment) = ntaxa(a)
Base.size(a::Alignment) = size(a.data)
Base.size(a::Alignment, d) = size(a.data, d)
Base.copy(a::Alignment) = Alignment(copy(a.names), copy(a.data))
Base.:(==)(a::Alignment, b::Alignment) = a.names == b.names && a.data == b.data

_decode(row) = String([_DEC[x+1] for x in row])
Base.getindex(a::Alignment, i::Integer) = _decode(view(a.data, i, :))
Base.getindex(a::Alignment, s::AbstractString) = a[_taxon(a, s)]
Base.iterate(a::Alignment, i=1) = i > ntaxa(a) ? nothing : (a.names[i] => a[i], i + 1)
Base.eltype(::Type{Alignment}) = Pair{String,String}

_taxon(a::Alignment, i::Integer) = i
function _taxon(a::Alignment, s::AbstractString)
    i = findfirst(==(s), a.names)
    i === nothing && throw(KeyError(s))
    return i
end
_taxa(a::Alignment, x::Colon) = 1:ntaxa(a)
_taxa(a::Alignment, x::AbstractVector{Bool}) = findall(x)
_taxa(a::Alignment, x) = [_taxon(a, i) for i in x]
_taxa(a::Alignment, x::Union{Integer,AbstractString}) = [_taxon(a, x)]

function Base.getindex(a::Alignment, t, s)
    rows = _taxa(a, t)
    Alignment(a.names[rows], a.data[rows, s isa Integer ? (s:s) : s])
end

"""
    sequences(aln) -> Vector{String}

The sequences of `aln` as strings (with gaps).
"""
sequences(a::Alignment) = [a[i] for i in 1:ntaxa(a)]

"""
    basefreqs(aln; ambiguous=true) -> Vector{Float64}

Empirical frequencies of A, C, G and T. Ambiguity codes are shared out
among the bases they allow, or ignored with `ambiguous = false` (as
`ape::base.freq` does); gaps and `N` are always ignored.
"""
function basefreqs(a::Alignment; ambiguous::Bool=true)
    f = zeros(4)
    for x in a.data
        m = x & 0x0f
        m == 0x0f && continue
        k = count_ones(m)
        k > 1 && !ambiguous && continue
        for b in 1:4
            (m >> (b - 1)) & 0x01 == 1 && (f[b] += 1 / k)
        end
    end
    s = sum(f)
    return s == 0 ? fill(0.25, 4) : f ./ s
end

# ------------------------------------------------------------------ site patterns

# Unique alignment columns and how often each occurs. Likelihood and
# parsimony are computed once per pattern; a bootstrap replicate is just a
# new set of weights.
struct _Patterns
    data::Matrix{UInt8}        # ntaxa × npatterns, gaps as 0x0f
    weight::Vector{Float64}
    site::Vector{Int}          # pattern of each original site
end

function _patterns(a::Alignment)
    n, m = size(a.data)
    seen = Dict{Vector{UInt8},Int}()
    cols = Vector{Vector{UInt8}}()
    weight = Float64[]
    site = Vector{Int}(undef, m)
    col = Vector{UInt8}(undef, n)
    for j in 1:m
        for i in 1:n
            col[i] = a.data[i, j] & 0x0f
        end
        k = get(seen, col, 0)
        if k == 0
            push!(cols, copy(col))
            push!(weight, 0.0)
            k = length(cols)
            seen[cols[k]] = k
        end
        weight[k] += 1
        site[j] = k
    end
    data = Matrix{UInt8}(undef, n, length(cols))
    for (k, c) in enumerate(cols)
        data[:, k] = c
    end
    return _Patterns(data, weight, site)
end

_npatterns(p::_Patterns) = size(p.data, 2)
_reweight(p::_Patterns, w::Vector{Float64}) = _Patterns(p.data, w, p.site)

# a bootstrap resample of the sites, as new pattern weights
function _resample(p::_Patterns, rng)
    w = zeros(length(p.weight))
    m = length(p.site)
    for _ in 1:m
        w[p.site[rand(rng, 1:m)]] += 1
    end
    return w
end

# constant (one base possible in every sequence), and parsimony-informative
# (at least two bases each seen unambiguously in at least two sequences)
function _siteclass(col::AbstractVector{UInt8})
    common = 0x0f
    cnt = (0, 0, 0, 0)
    for x in col
        m = x & 0x0f
        common &= m
        if _isbase(m)
            b = _index(m)
            cnt = Base.setindex(cnt, cnt[b] + 1, b)
        end
    end
    return common != 0, count(>=(2), cnt) >= 2
end

# ------------------------------------------------------------------ reading

"""
    readalignment(path) -> Alignment

Read aligned sequences from a FASTA, PHYLIP (sequential or interleaved) or
Nexus (`DATA`/`CHARACTERS` block) file; the format is detected from the
contents.
"""
function readalignment(path::AbstractString)
    s = read(path, String)
    st = lstrip(s)
    startswith(st, '>') && return _parse_fasta(s)
    startswith(uppercase(st), "#NEXUS") && return _parse_nexus_matrix(s)
    occursin(r"^\s*\d+\s+\d+", st) && return _parse_phylip(s)
    throw(ArgumentError("$path: not FASTA, PHYLIP or Nexus"))
end

"""    readfasta(path) -> Alignment  read an aligned FASTA file."""
readfasta(path::AbstractString) = _parse_fasta(read(path, String))
"""    readphylip(path) -> Alignment  read a PHYLIP file (sequential or interleaved)."""
readphylip(path::AbstractString) = _parse_phylip(read(path, String))

function _parse_fasta(s::AbstractString)
    names = String[]
    seqs = String[]
    buf = IOBuffer()
    for line in eachline(IOBuffer(s))
        line = rstrip(line)
        if startswith(line, '>')
            isempty(names) || push!(seqs, String(take!(buf)))
            push!(names, String(strip(line[2:end])))
        elseif !isempty(names)
            print(buf, filter(!isspace, line))
        end
    end
    isempty(names) || push!(seqs, String(take!(buf)))
    # keep the first word as the name when that keeps names unique (FASTA
    # headers often carry a description after the id)
    ids = [isempty(n) ? n : first(split(n)) for n in names]
    return Alignment(allunique(ids) ? ids : names, seqs)
end

function _parse_phylip(s::AbstractString)
    lines = filter(!isempty, strip.(split(s, '\n')))
    hdr = split(lines[1])
    n, m = parse(Int, hdr[1]), parse(Int, hdr[2])
    r = _phylip_sequential(lines[2:end], n, m)
    r === nothing && (r = _phylip_interleaved(lines[2:end], n, m))
    r === nothing && throw(ArgumentError("PHYLIP: could not read $n sequences of $m sites"))
    return Alignment(r...)
end

# "name  ACGT…" (relaxed: any whitespace after the name) or, when there is
# no whitespace, the strict format's 10-character name field
function _phylip_nameline(line)
    parts = split(line; limit=2)
    length(parts) == 2 && return String(parts[1]), filter(!isspace, parts[2])
    return String(strip(line[1:min(10, end)])), filter(!isspace, line[min(11, end + 1):end])
end

function _phylip_sequential(lines, n, m)
    names, seqs = String[], String[]
    i = 1
    for _ in 1:n
        i > length(lines) && return nothing
        name, q = _phylip_nameline(lines[i]); i += 1
        while length(q) < m && i <= length(lines)
            q *= filter(!isspace, lines[i]); i += 1
        end
        length(q) == m || return nothing
        push!(names, name); push!(seqs, q)
    end
    return names, seqs
end

function _phylip_interleaved(lines, n, m)
    length(lines) < n && return nothing
    names, seqs = String[], String[]
    for k in 1:n
        name, q = _phylip_nameline(lines[k])
        push!(names, name); push!(seqs, q)
    end
    for (j, line) in enumerate(lines[n+1:end])
        k = mod1(j, n)
        seqs[k] *= filter(!isspace, line)
    end
    all(q -> length(q) == m, seqs) || return nothing
    return names, seqs
end

function _parse_nexus_matrix(s::AbstractString)
    m = match(r"\bmatrix\b(.*?);"is, s)
    m === nothing && throw(ArgumentError("Nexus file has no MATRIX"))
    order = String[]
    seqs = Dict{String,IOBuffer}()
    for line in split(m.captures[1], '\n')
        line = strip(replace(line, r"\[.*?\]" => ""))
        isempty(line) && continue
        if startswith(line, '\'')
            j = findnext('\'', line, 2)
            name, rest = line[2:j-1], line[j+1:end]
        else
            parts = split(line; limit=2)
            length(parts) < 2 && continue
            name, rest = parts
        end
        name = String(name)
        if !haskey(seqs, name)
            push!(order, name)
            seqs[name] = IOBuffer()
        end
        print(seqs[name], filter(!isspace, rest))
    end
    return Alignment(order, [String(take!(seqs[n])) for n in order])
end

# ------------------------------------------------------------------ writing

"""
    writefasta(path, aln; width=60)
    writefasta(io, aln; width=60)

Write `aln` as FASTA, wrapping sequences at `width` characters (`0`: no
wrapping).
"""
function writefasta(io::IO, a::Alignment; width::Integer=60)
    for i in 1:ntaxa(a)
        println(io, '>', a.names[i])
        s = a[i]
        if width <= 0
            println(io, s)
        else
            for k in 1:width:length(s)
                println(io, s[k:min(end, k + width - 1)])
            end
        end
    end
end
writefasta(path::AbstractString, a::Alignment; kwargs...) = (open(io -> writefasta(io, a; kwargs...), path, "w"); path)

"""
    writephylip(path, aln)
    writephylip(io, aln)

Write `aln` as relaxed sequential PHYLIP (names of any length, followed by
a space).
"""
function writephylip(io::IO, a::Alignment)
    println(io, ntaxa(a), ' ', nsites(a))
    w = maximum(length, a.names; init=0)
    for i in 1:ntaxa(a)
        println(io, rpad(replace(a.names[i], ' ' => '_'), w), "  ", a[i])
    end
end
writephylip(path::AbstractString, a::Alignment) = (open(io -> writephylip(io, a), path, "w"); path)
