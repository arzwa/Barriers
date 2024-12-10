
@kwdef struct MainlandIslandModel{T,V,A<:Architecture}
    arch :: A
    m    :: T
    p̄    :: Vector{T} = zeros(length(arch))
    N    :: V   # may be a vector, if we want to include Nₑ heterogeneity
end

struct Equilibrium{M<:MainlandIslandModel,T}
    model :: M
    Ep    :: Vector{T} 
    Epq   :: Vector{T}
end

function Equilibrium(model::MainlandIslandModel; 
        p0=ones(length(model.arch)), tol=1e-9)
    @unpack arch, p̄, N, m = model
    @unpack R, loci = arch
    pq0 = p0 .* (1 .- p0)
    Ep, Epq = fixed_point_iteration(m, loci, R, N, p̄, p0, pq0, tol=tol)
    Equilibrium(model, Ep, Epq)
end

getNe(N::Real, _) = N
getNe(N, i) = N[i]

function fixed_point_iteration(m, loci, R, N, p̄, p, pq; tol=1e-9)
    L = length(loci) 
    while true
        gs  = gffs(m, loci, p̄, R, p, pq) 
        xs  = map(i->predict(loci[i], getNe(N,i), p̄[i], m*gs[i]), 1:L)
        Ep  = first.(xs)
        Epq = last.(xs)
        norm(Ep .- p) < tol && return Ep, Epq
        p = Ep; pq = Epq
    end
end

gffs(m, loci, p̄, R, p, pq) = map(i->gff(m, loci, p̄, R, p, pq, i), 1:length(loci))

function gff(m, loci, p̄, R, p, pq, i)
    lg = 0.0
    for j=1:length(loci)
        i == j && continue  # it does matter! We overpredict the beneficial
        # allele frequency when we include the focal locus' own effect in the
        # mₑ calculation
        lg += locuseffect(loci[j], p[j], pq[j], m, p̄[j], R[i,j])
    end
    exp(-lg)
end

function locuseffect(locus, p, pq, m, p̄, r)
    @unpack s, h = locus
    N = s*h*(p-p̄) - s*(1-2h)*(p̄*(1-p) - pq)
    D = m + r + s*(h - (1-p) + 2*(1-2h)*pq)
    return N/D
end

function predict(locus, Nₑ, p̄, mₑ)
    @unpack s, h, u = locus
    # p is locally beneficial, q is locally deleterious
    d = Wright(-Nₑ*s, Nₑ*(mₑ*(p̄ + u)), Nₑ*(mₑ*(1-p̄) + u), h) 
    # mean(d) is the frequency of the locally beneficial allele
    mean(d), expectedpq(d)
end

function gff(eqmodel::Equilibrium, x)
    @unpack model, Ep, Epq = eqmodel
    @unpack m, p̄, arch = model
    @unpack xs, loci = arch
    lg = 0.0
    for j=1:length(loci)
        lg += locuseffect(loci[j], Ep[j], Epq[j], m, p̄[j], recrate(abs(x - xs[j])))
    end
    exp(-lg)
end

me(M::Equilibrium, x) = M.model.m*gff(M, x)

