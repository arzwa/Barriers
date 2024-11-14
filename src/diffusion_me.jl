"""
    MILocus

- `s`  selection coefficient *against island-locally del. allele*
- `h`  dominance coefficient
- `u`  mutation rate (symmetric mutation assumed, could be relaxed)
- `p̄`  mainland allele frequency for the island-locally del. allele
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
        xs  = map(i->predict(loci[i], m*gs[i]), 1:length(loci))
        Ep  = first.(xs)
        Epq = last.(xs)
        norm(Ep .- p) < tol && return Ep, Epq
        p = Ep; pq = Epq
    end
end

function predict(locus::MILocus, mₑ)
    @unpack Ne, s, h, u, p̄ = locus
    # p is locally beneficial, q is locally deleterious
    d = Wright(-Ne*s, Ne*(mₑ*(p̄ + u)), Ne*(mₑ*(1-p̄) + u), h) 
    # mean(d) is the frequency of the locally beneficial allele
    mean(d), expectedpq(d)
end

gffs(m, loci, R, p, pq) = map(i->gff(m, loci, R, p, pq, i), 1:length(loci))

function gff(m, loci, R, p, pq, i)
    lg = 0.0
    for j=1:length(loci)
        i == j && continue  # it does matter! We overpredict the beneficial
        # allele frequency when we include the focal locus' own effect in the
        # mₑ calculation
        lg += locuseffect(loci[j], p[j], pq[j], m, R[i,j])
    end
    exp(-lg)
end

function locuseffect(locus, p, pq, m, r)
    @unpack s, h, p̄ = locus
    N = s*h*(p-p̄) - s*(1-2h)*(p̄*(1-p) - pq)
    D = m + r + s*(h - (1-p) + 2*(1-2h)*pq)
    return N/D
end

function gff(model::MIModel, x)
    @unpack m, xs, loci, Ep, Epq = model
    lg = 0.0
    for j=1:length(loci)
        lg += locuseffect(loci[j], Ep[j], Epq[j], m, recrate(x, xs[j]))
    end
    exp(-lg)
end

me(model, x) = model.m*gff(model, x)

