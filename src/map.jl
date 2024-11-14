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

maplength(m::GeneticMap) = last(m.gen2phys.itp.knots[1])
physlength(m::GeneticMap) = last(m.phys2gen.itp.knots[1])

distance(d::GeneticMap, x1, x2) = abs(d[x1] - d[x2])
recrate(d::GeneticMap, x1, x2) = recrate(distance(d, x1, x2))
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
    
"""
    WindowedChromosome(genetic_map::GeneticMap, winsize, winstep)
"""
struct WindowedChromosome{T,V,U,W}
    data::Vector{Tuple{UnitRange{T},Tuple{V,V}}}
    geneticmap::GeneticMap{U,W}
end

function Base.show(io::IO, d::WindowedChromosome)
    x = (length(d), physlength(d)/Mb, maplength(d)*100)
    write(io, @sprintf("WindowedChromosome(#%d, %.1fMb, %.1fcM)", x...))
end

"""
    WindowedChromosome(gmap, winsize, winstep)

Construct a WindowedChromosome using regular windowsize stepsize.
"""
function WindowedChromosome(gmap::GeneticMap, winsize::Int, winstep::Int) 
    @assert winsize % winstep == 0 "`winsize` should be a multiple of `winstep`"
    chr_end = physlength(gmap)
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
    WindowedChromosome(windows, gmap) 
end

"""
    WindowedChromosome(gmap, windows::Vector)

Construct a WindowedChromosome using given windows (can be irregular).
"""
function WindowedChromosome(gmap::GeneticMap, windows::Vector)
    xs = map(windows) do win
        x1, x2 = extrema(win)
        (x1:x2, (gmap[x1], gmap[x2]))
    end
    WindowedChromosome(xs, gmap)
end

Base.getindex(d::WindowedChromosome, i) = d.data[i]
Base.length(d::WindowedChromosome) = length(d.data)
Base.lastindex(d::WindowedChromosome) = length(d)

maplength(d::WindowedChromosome) = maplength(d.geneticmap)
physlength(d::WindowedChromosome) = physlength(d.geneticmap)

# Coarse recombination rates between windows (recombination rates between
# window midpoints).
function window_recrates(d::WindowedChromosome)
    n = length(d)
    R = zeros(n, n)
    for i=2:n
        x0, x1 = d[i][2]
        xmid = (x0 + x1)/2
        for j=1:i-1
            x0, x1 = d[j][2]
            R[i,j] = R[j,i] = recrate(abs.((x0 + x1)/2 - xmid))
        end
    end
    return R
end

"""
    me_profile(model, d::WindowedChromosome)

Calculate expected `mₑ` in windows under the given model.
"""
function me_profile(model, d::WindowedChromosome)
    map(1:length(d)) do i
        (window, (x1, x2)) = d[i]
        Eme, _ = quadgk(x->me(model, x), x1, x2)
        window => Eme/(x2 - x1)
    end 
end

function me_profile_vec(model, d::WindowedChromosome)
    map(1:length(d)) do i
        (_, (x1, x2)) = d[i]
        Eme, _ = quadgk(x->me(model, x), x1, x2)
        Eme/(x2 - x1)
    end 
end

function me_profile_vec_local!(
        model, d::WindowedChromosome, 
        mes, loci, j)
    @unpack n = model
    # recompute focal window me
    (_, (x1, x2)) = d[j]
    Eme, _ = quadgk(x->me(model, x), x1, x2)
    mes[j] = Eme/(x2 - x1)
    # examine leftward windows
    left = 0; i = 1
    while left < n && j - i > 0
        (_, (x1, x2)) = d[j-i]
        Eme, _ = quadgk(x->me(model, x), x1, x2)
        mes[j-i] = Eme/(x2 - x1) 
        left += length(loci[j-i])
        i += 1
    end
    # examine rightward windows
    rght = 0; i = 1
    while rght < n && j + i < length(d)
        (_, (x1, x2)) = d[j+i]
        Eme, _ = quadgk(x->me(model, x), x1, x2)
        mes[j+i] = Eme/(x2 - x1) 
        left += length(loci[j+i])
        i += 1
    end
    return mes
end

"""
    plotcoords

Takes output from `me_profile` and returns plotting coordinates.
"""
function plotcoords(zs)
    xs = vcat(map(x->[extrema(x)...], first.(zs))...)
    ys = vcat(map(x->[x, x], last.(zs))...)
    xs, ys
end


