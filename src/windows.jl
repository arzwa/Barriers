"""
    GeneticMap(x::Vector{Tuple{<:Int,<:Real}})

This interfaces between a physical map in units of **basepairs** and a genetic
map in units of **Morgans** (i.e. expected number of crossovers). This is for a
single DNA molecule (chromosome). To construct an instance, provide a list with
tuples of physical locations and corresponding map positions, including the
origin of the map `(0, 0.0)`.

When indexing with an integer `i`, the genetic distance from the start of the
chromosome to base `i` is returned (in Morgans). When indexing with a real
number `x`, the closest basepair to the map position `x`M is returned. 

```julia
julia> gmap = GeneticMap([(0,0.0), (200000, 0.002), (1000000, 0.004)])
```
"""
struct GeneticMap{T,V}
    phys2gen::T
    gen2phys::V
end

Base.getindex(m::GeneticMap, x::Real) = round(Int64, m.gen2phys[x])
Base.getindex(m::GeneticMap, x::Int) = m.phys2gen[x]
   
function GeneticMap(genetic_map)
    phys2gen = linear_interpolation(first.(genetic_map), last.(genetic_map))
    gen2phys = linear_interpolation(last.(genetic_map), first.(genetic_map))
    GeneticMap(phys2gen, gen2phys)
end

distance(d::GeneticMap, x1, x2) = abs(d[x1] - d[x2])
recrate(d::GeneticMap, x1, x2) = haldanemapfun(distance(d, x1, x2))
recrate(x1, x2) = recrate(abs(x1 - x2))
recrate(d) = 0.5*(1-exp(-2d))

function recrates(xs::Vector{T}) where T
    L = length(xs)
    R = zeros(L, L)
    for i=1:L
        for j=1:i-1
            R[i,j] = R[j,i] = recrate(xs[i], xs[j])
        end
    end
    return R
end
    
"""
    WindowedChromosome(genetic_map::Vector{Tuple}, winsize)

`genetic_map` should be a list of pairs (physical position [b], genetic
position [M]). `winsize` is also in b. Note that we expect (0,0.0) to be in the
map.
"""
struct WindowedChromosome{T,V,U,W}
    data::Vector{Tuple{UnitRange{T},Tuple{V,V}}}
    winsize::Int64
    winstep::Int64
    nonoverlapping::StepRange{Int64,Int64}
    geneticmap::GeneticMap{U,W}
end

function Base.show(io::IO, d::WindowedChromosome)
    x = (length(d), d.winsize/kb, d.winstep/kb, physlength(d)/Mb, maplength(d)*100)
    write(io, @sprintf("WindowedChromosome(#%d, %.1fkb, %.1fkb, %.1fMb, %.1fcM)", x...))
end

function WindowedChromosome(genetic_map, winsize, winstep) 
    @assert issorted(genetic_map)
    @assert winsize % winstep == 0
    gmap = GeneticMap(genetic_map)
    chr_end = last(genetic_map)[1]
    x1, x2 = 1, winsize
    y1, y2 = gmap[x1], gmap[x2]
    windows = [(x1:x2, (y1,y2))]
    while true
        x1 += winstep 
        x2 = min(chr_end, x2 + winstep)
        y1 = gmap[x1]
        y2 = gmap[x2]
        push!(windows, (x1:x2, (y1,y2)))
        x2 >= chr_end && break
    end
    Δ = winsize÷winstep
    nonoverlapping = 1:Δ:length(windows)
    WindowedChromosome(windows, winsize, winstep, nonoverlapping, gmap) 
end

Base.getindex(d::WindowedChromosome, i) = d.data[i]
Base.length(d::WindowedChromosome) = length(d.data)
Base.lastindex(d::WindowedChromosome) = length(d)

maplength(d::WindowedChromosome) = d[end][2][2]
physlength(d::WindowedChromosome) = d[end][1][end]

