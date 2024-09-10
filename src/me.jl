# We want to calculate an expected m_e in windows.
# What we should get as input is the selected architecture: the selected loci
# and their map positions, the rest should be done under the hood.
struct AeschbacherModel{T}
    m::T          # migration rate
    s::Vector{T}  # selection coefficients
    xs::Vector{T}  # map positions of the selected loci
    function AeschbacherModel(m::T, s::Vector{T}, xs::Vector{T}) where T 
        @assert issorted(xs)
        new{T}(m, s, xs)
    end
end

# calculate me at **map** position x
function me(model::AeschbacherModel, x)
    @unpack m, s, xs = model
    loggff = 0.0
    function _recurse1(a, i)
        (i == length(xs) + 1 || xs[i] > x) && return 0.0
        b = _recurse1(a, i+1)
        r = recrate(xs[i], x)
        loggff -= log(1 - s[i]/(r + b))
        return a - b  # Aeschbacher expressions are for positive s...
    end
    function _recurse2(a, i)
        (i == 0 || xs[i] < x) && return 0.0
        b = _recurse2(a, i-1)
        r = recrate(xs[i], x)
        loggff -= log(1 - s[i]/(r + b))
        return a - b
    end
    _recurse1(0.0, 1)
    _recurse2(0.0, length(s))
    return m*exp(loggff)
end

function me_profile(model, d::WindowedChromosome)
    map(1:length(d)) do i
        (window, (x1, x2)) = d[i]
        Eme, _ = quadgk(x->me(model, x), x1, x2)
        window => Eme/(x2 - x1)
    end
end

struct MILocus{T}
    s  :: T
    h  :: T
    u  :: T
    p̄  :: T
    Ne :: T
end

"""
    MainlandIslandModel

Polygenic barrier model as in Zwaenepoel et al. (2024) and Sachdeva (2022).
"""
struct MIModel{T}
    m    :: T
    xs   :: Vector{T}
    loci :: Vector{MILocus{T}}
    Ep   :: Vector{T}
    Epq  :: Vector{T}
end

function MIModel(m::T, xs::Vector{T}, loci; p=ones(length(xs)), tol=1e-9) where T
    @assert length(xs) == length(loci)
    R = recrates(xs)
    Ep, Epq = fixed_point_iteration(m, loci, R, p, p .* (1 .- p), tol=tol)
    MIModel(m, xs, loci, Ep, Epq)
end

function fixed_point_iteration(m, loci, R, p, pq; tol=1e-9)
    while true
        gs  = gffs(m, loci, R, p, pq) 
        @info gs
        xs  = map(i->predict(loci[i], m*gs[i]), 1:length(loci))
        Ep  = first.(xs)
        Epq = last.(xs)
        norm(Ep .- p) < tol && return Ep, Epq
        p = Ep; pq = Epq
    end
end

function predict(locus::MILocus, mₑ)
    @unpack Ne, s, h, u, p̄ = locus
    d = Wright(Ne*s, Ne*(mₑ*(1-p̄) + u), Ne*(mₑ*p̄ + u), h) 
    mean(d), expectedpq(d)
end

gffs(m, loci, R, p, pq) = map(i->gff(m, loci, R, p, pq, i), 1:length(loci))

function gff(m, loci, R, p, pq, i)
    lg = 0.0
    for j=1:length(loci)
        lg += locuseffect(loci[j], p[j], pq[j], m, R[i,j])
    end
    exp(-lg)
end

function locuseffect(locus, p, pq, m, r)
    @unpack s, h, p̄ = locus
    N = -s*h*(p-p̄) + s*(1-2h)*(p̄*(1-p) - pq)
    D = m + r + -s*(h - (1-p) + 2*(1-2h)*pq)
    return N/D
end

function me(model::MIModel, x)
    @unpack m, xs, loci, Ep, Epq = model
    lg = 0.0
    for j=1:length(loci)
        lg += locuseffect(loci[j], Ep[j], Epq[j], m, recrate(x, xs[j]))
    end
    m*exp(-lg)
end
