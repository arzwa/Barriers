@kwdef struct MainlandIslandModel{T,V,A<:Architecture}
    arch :: A
    m    :: T
    p̄    :: Vector{T} = zeros(length(arch))
    N    :: V   # may be a vector, if we want to include Nₑ heterogeneity
end

struct Equilibrium{M,T}
    model :: M
    Ep    :: T
    Epq   :: T
end

function Equilibrium(model::MainlandIslandModel; 
        p0=ones(length(model.arch)), tol=1e-9)
    @unpack arch, p̄, N, m = model
    @unpack R, loci = arch
    pq0 = p0 .* (1 .- p0)
    Ep, Epq = fixed_point_iteration(m, loci, R, N, p̄, p0, pq0, tol=tol)
    Equilibrium(model, Ep, Epq)
end

function distributions(eq::Equilibrium)
    @unpack arch, p̄, N, m = eq.model
    @unpack R, loci = arch
    gs = gffs(m, loci, p̄, R, eq.Ep, eq.Epq) 
    map(1:length(loci)) do i
        @unpack s, h, u = loci[i]
        d = Wright(-N*s, N*(m*gs[i]*p̄[i] + u), N*(m*gs[i]*(1-p̄[i]) + u), h) 
    end
end

getNe(N::Real, _) = N
#getNe(N, i) = N[i]

function fixed_point_iteration(m, loci, R, N, p̄, p, pq; tol=1e-9)
    L = length(loci) 
    while true
        gs  = gffs(m, loci, p̄, R, p, pq) 
        xs  = map(i->predict(loci[i], getNe(N,i), p̄[i], m*gs[i]), 1:L)
        Ep  = first.(xs)
        Epq = last.(xs)
        norm(Ep .- p) < tol && return (Ep, Epq)
        p  = Ep
        pq = Epq
    end
end

gffs(m, loci, p̄, R, p, pq) = map(
    i->gff(m, loci, p̄, R, p, pq, i), 1:length(loci))

function gff(m, loci, p̄, R, p, pq, i)
    lg = 0.0
    for j=1:length(loci)
        i == j && continue  # it does matter! We overpredict the beneficial
        # allele frequency when we include the focal locus' own effect in the
        # mₑ calculation
        lg += locuseffect(loci[j], p[j], pq[j], m, p̄[j], R[i,j])
    end
    exp(lg)
end

function locuseffect(locus, p, pq, m, p̄, r)
    @unpack s, h = locus
    q = 1 - p
    q̄ = 1 - p̄
    #sa = -s*h
    #sb = -s*(1 - 2h)
    #q  = 1 - p
    #N  = sa*(p - p̄) + sb*(pq - p̄*q)
    #D  = m + r - sa*(p - q) - sb*(2pq - q)  
    # below is wrong? NOT the same result when dominance??
    N = -s*h*(q̄ - q) + s*(1-2h)*(p̄*q - pq)
    D = m + r + s*(h - (1-p) + 2*(1-2h)*pq)
    # haploid -> sₑ = 2s, h=1/2 -> m + r + 2s(1/2 - (1-p)) = m + r - s(1 - 2p)
    return N/D
end

function predict(locus, Nₑ, p̄, mₑ)
    @unpack s, h, u = locus
    # p̄ = 0
    # p is locally beneficial, q is locally deleterious
    # migration introduces q
    q̄ = 1 - p̄
    d = Wright(-Nₑ*s, Nₑ*(mₑ*p̄ + u), Nₑ*(mₑ*q̄ + u), h) 
    # @info -Nₑ*s, Nₑ*(mₑ*p̄ + u), Nₑ*(mₑ*q̄ + u), mean(d)
    # mean(d) is the frequency of the locally beneficial allele
    mean(d), expectedpq(d)
end

function gff(eqmodel::Equilibrium, x)
    @unpack model, Ep, Epq = eqmodel
    @unpack m, p̄, arch = model
    @unpack xs, loci = arch
    lg = 0.0
    for j=1:length(loci)
        rij = recrate(abs(x-xs[j]))
        le = locuseffect(loci[j], Ep[j], Epq[j], m, p̄[j], rij)
        lg += le
    end
    exp(lg)
end

me(M::Equilibrium, x) = M.model.m*gff(M, x)


# Unlinked haploid. mₑ = m*exp(-2Ls𝔼p). Assume mainland fixed (we predict
# divergence).
struct UnlinkedHaploidMainlandIsland{T}
    L :: T
    s :: T
    m :: T
    u :: T
    N :: T
end
UnlinkedHaploidMainlandIsland(args...) = 
    UnlinkedHaploidMainlandIsland(promote(args...)...)

function Equilibrium(model::UnlinkedHaploidMainlandIsland; p0=1.0, tol=1e-9)
    @unpack s, m, u, N, L = model
    r̄ = 0.5
    mₑ = m
    while true
        # s is effective selection! i.e. haploid selection sₕ => s = 2sₕ
        Ep, _ = predict((s=s, u=u, h=0.5), N, 0.0, mₑ)
        _mₑ = m*exp(-L*s*(1/2)*Ep/r̄)  # s here is twice the haploid selection coeff! 
        # we need exp(-2Lsₑhₑ𝔼p) = exp(-Lsₑ𝔼p) 
        abs(mₑ - _mₑ) < tol && return Equilibrium(model, Ep, Ep*(1-Ep))
        mₑ = _mₑ
    end
end

function UnlinkedEquilibriumApproxHom(model::MainlandIslandModel; tol=1e-5)
    # assumes all loci are identical (`Hom` for homogeneous)
    @unpack arch, m, N = model
    L = length(arch)
    M_ = UnlinkedHaploidMainlandIsland(L, arch[1].s, m, arch[1].u, N)
    eq = Equilibrium(M_; tol=tol)
    Equilibrium(model, fill(eq.Ep, L), fill(eq.Epq, L))
end

function ApproxEquilibrium(model::MainlandIslandModel; tol=1e-5)
    @unpack arch, m, N = model
    L = length(arch)
    s = arch[1].s
    u = arch[1].u
    mₑ = m
    while true
        # s is effective selection! i.e. haploid selection sₕ => s = 2sₕ
        Ep, _ = predict((s=s, u=u, h=0.5), N, 0.0, mₑ)
        #r̄ = _denom_approx(arch, m, s, Ep)
        r̄ = 0.5
        _mₑ = m*exp(-(L-1)*s*(1/2)*Ep/r̄)  
        # s here is twice the haploid selection coeff! 
        # we need exp(-2Lsₑhₑ𝔼p) = exp(-Ls𝔼p) 
        abs(mₑ-_mₑ) < tol && return Equilibrium(model, fill(Ep, L), fill(Ep*(1-Ep),L))
        mₑ = _mₑ
    end
end

# crude way to account somewhat for linkage
function _denom_approx(arch, m, s, Ep)
    L = length(arch)
    i = L÷2
    D = 0.0
    for j=1:i-1
        D += 1/(m + arch.R[i,j] + s*(2Ep - 1))
    end
    for j=i+1:L
        D += 1/(m + arch.R[i,j] + s*(2Ep - 1))
    end
    return (L-1)/D
end



