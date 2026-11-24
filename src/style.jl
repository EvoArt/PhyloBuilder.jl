# Terminal styling shared by the object displays (show.jl) and the progress
# bars (progress.jl), in the same look as DenoiseDNA.jl and FindRNA.jl.
#
# Colours are xterm-256 indices: they look the same in every modern terminal,
# whereas the 16 named colours are re-themed per terminal. Escape codes are
# written only when the IO has :color, so sprint(), logs and files get plain
# text.

const _PAL = (head=117, key=245, rule=240, dim=242, good=114, warn=221, amber=214, bad=203,
              spark=117, tail=239, hi=255,
              A=114, C=75, G=214, T=203, N=245, gap=239)

_hascolor(io) = get(io, :color, false)::Bool

# `s` in an xterm-256 colour (and/or bold), or plain `s` when the IO has no colour
function _sty(io, s; color=nothing, bold::Bool=false)
    _hascolor(io) || return string(s)
    pre = (bold ? "\e[1m" : "") * (color === nothing ? "" : "\e[38;5;$(color)m")
    return isempty(pre) ? string(s) : string(pre, s, "\e[0m")
end

# 1234567 -> "1,234,567"
function _fmt(n::Integer)
    s = string(abs(n))
    k = length(s)
    b = IOBuffer()
    n < 0 && print(b, '-')
    for (i, c) in enumerate(s)
        print(b, c)
        r = k - i
        r > 0 && r % 3 == 0 && print(b, ',')
    end
    return String(take!(b))
end
_fmt(x::Real) = isinteger(x) && abs(x) < 1e15 ? _fmt(Int(x)) : _g(x)

_g(x) = isnan(x) ? "NaN" : isinf(x) ? (x > 0 ? "Inf" : "-Inf") :
        x == 0 ? "0" : abs(x) >= 1e4 || abs(x) < 1e-2 ? @sprintf("%.2e", x) : @sprintf("%.3g", x)
