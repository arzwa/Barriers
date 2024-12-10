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

function recrate(d::WindowedChromosome, win)
    xx, (d0, d1) = d[win]
    r = recrate(d.geneticmap, xx[1], xx[2])
    r, length(xx)
end


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


