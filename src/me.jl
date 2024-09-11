"""
    AeschbacherModel(m, s::Vector, xs::Vector)

- `m`  migration rate
- `s`  vector of selection coefficients (haploid selection)
- `xs` vector of map positions (in Morgan)
"""
struct AeschbacherModel{T}
    m :: T          # migration rate
    s :: Vector{T}  # selection coefficients
    xs:: Vector{T}  # map positions of the selected loci
    n :: Int64      # consider n MSPs on either side
    function AeschbacherModel(m::T, s::Vector{T}, xs::Vector{T}; n=100) where T 
        new{T}(m, s, sort(xs), n)
    end
end

# calculate me at **map** position x
function me(model::AeschbacherModel, x)
    @unpack m, s, xs = model
    loggff = 0.0
    function _recursel(a, i)
        (i == length(xs) + 1 || xs[i] > x) && return 0.0
        b = _recursel(a, i+1)
        r = recrate(xs[i], x)
        loggff -= log(1 - s[i]/(r + b))
        return a - b  # Aeschbacher expressions are for positive s...
    end
    function _recurser(a, i)
        (i == 0 || xs[i] < x) && return 0.0
        b = _recurser(a, i-1)
        r = recrate(xs[i], x)
        loggff -= log(1 - s[i]/(r + b))
        return a - b
    end
    left = findlast(z->z < x, xs)
    left = isnothing(left) ? 1 : left
    xl = max(1, left - model.n + 1)
    xr = min(length(s), left + model.n)
    _recursel(0.0, xl)
    _recurser(0.0, xr)
    return m*exp(loggff)
end

"""
    MILocus

- `s`  selection coefficient
- `h`  dominance coefficient
- `u`  mutation rate (symmetric mutation assumed, could be relaxed)
- `p̄`  mainland allele frequency for the island-locally beneficial allele
- `Ne` effective population size
"""
struct MILocus{T}
    s  :: T
    h  :: T
    u  :: T
    p̄  :: T
    Ne :: T
end

"""
   MIModel(m, xs, loci::Vector{MILocus}; p, tol)

Polygenic barrier model as in Zwaenepoel et al. (2024) and Sachdeva (2022).

- `m`    migration rate
- `xs`   map positions of selected loci
- `loci` vector of selected loci
- `p`    initial allele frequencies for fixed point iteration (default
         `ones(L)`, i.e. secondary contact scenario).
- `tol`  tolerance for the fixed point iteration (default `1e-9`)
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

