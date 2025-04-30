"""
    GeneticMap(x::Vector{Tuple{<:Int,<:Real}})

This interfaces between a physical map in units of **basepairs** and a genetic
map in units of **Morgans** (i.e. expected number of crossovers). This is for a
single DNA molecule (chromosome). To construct an instance, provide a list with
tuples of physical locations and corresponding map positions, including the
origin of the map `(0, 0.0)`.

When indexing with an integer `i`, the genetic distance from the start of the
chromosome to base `i` is returned (in Morgans). When indexing with a real
number `x`, the closest basepair to the map position `x`M is returned. We
interpolate linearly between markers.

```julia
julia> gmap = GeneticMap([(0,0.0), (200000, 0.002), (1000000, 0.004)])
```
"""
struct GeneticMap{T,V}
    data::Vector{Tuple{Int64,Float64}}
    phys2gen::T
    gen2phys::V
end

function GeneticMap(genetic_map)
    @assert issorted(genetic_map) "Map positions should be sorted"
    if !(genetic_map[1] == (1, 0.0))
        genetic_map = [(1, 0.0) ; genetic_map]
    end
    phys2gen = linear_interpolation(first.(genetic_map), last.(genetic_map))
    gen2phys = linear_interpolation(last.(genetic_map), first.(genetic_map))
    GeneticMap(genetic_map, phys2gen, gen2phys)
end

linearmap(physlen, maplen) = GeneticMap([(1,0.0), (physlen,maplen)])

Base.getindex(m::GeneticMap, x::Real) = round(Int64, m.gen2phys[x])
Base.getindex(m::GeneticMap, x::Int) = m.phys2gen[x]
Base.show(io::IO, m::GeneticMap) = write(io, 
    @sprintf("GeneticMap(%.1fMb, %.1fcM)", physlength(m)/Mb, maplength(m)*100))

maplength( m::GeneticMap) = last(m.gen2phys.itp.knots[1])
physlength(m::GeneticMap) = last(m.phys2gen.itp.knots[1])

recrate( d::GeneticMap, x1::Int, x2::Int) = recrate(distance(d, x1, x2))
distance(d::GeneticMap, x1::Int, x2::Int) = abs(d[x1] - d[x2])
recrate(d)  = 0.5*(1-exp(-2d))
distance(r) = -0.5log(1-2r)
    
"""
    recrates(xs::Vector)
    recrates(g::GeneticMap, xs::Vector)

Calculate recombination rate matrix between the map positions `xs`.
"""
function recrates(xs::AbstractVector{T}) where T
    L = length(xs)
    R = zeros(L, L)
    for i=1:L
        for j=1:i-1
            R[i,j] = R[j,i] = recrate(abs(xs[i] - xs[j]))
        end
    end
    return R
end

recrates(g::GeneticMap, xs::Vector{Int}) = recrates([g[x] for x in xs])

function average_rr(g::GeneticMap; scale::Int64=1kb)
    xs = 1:scale:physlength(g)
    rr = 0.0 
    for i=2:length(xs)
        rr += distance(g, xs[i], xs[i-1])/(xs[i] - xs[i-1])
    end
    rr / (length(xs)-1)  # M/bp
end

"""
    extrapolate_to(m::GeneticMap, chr_end::Int64; scale=10kb)

Extrapolate the genetic map to a physical position `chr_end`, using the avreage
recombination rate estimated over a scale of `scale`.
"""
function extrapolate_to(m::GeneticMap, chr_end::Int64; scale=10kb)
    @assert physlength(m) < chr_end "Map is longer than provided new endpoint"
    rr = average_rr(m, scale=scale)
    ml = chr_end * rr  # maplength
    GeneticMap([m.data; (chr_end, ml)])
end

