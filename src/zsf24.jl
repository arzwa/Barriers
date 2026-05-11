
@with_kw struct ZSF24{T,V,R<:RecombinationMap,N}
    m  :: T
    s  :: Vector{T}
    h  :: Vector{T}
    xs :: Vector{V}
    rm :: R
    Ne :: N
    u  :: T
    p̄  :: Vector{T}
    mode = 1
end

function gff(model::ZSF24, Ep, Epq)
    @unpack s, h, m, rm, xs, p̄, mode = model
    L = length(xs)
    sp = s .* Ep
    map(1:L) do k
        gff(rm, m, xs, s, h, Ep, Epq, p̄, k, mode)
    end
end

function gff(rm, m, xs, s, h, p, pq, p̄, k, mode=1)
    loggff = 0.0
    for j=1:length(xs)
        j == k && continue  # it does matter! We overpredict the beneficial
        rjk = recrate(rm, xs[j], xs[k])
        loggff += locuseffect(s[j], h[j], p[j], pq[j], m, p̄[j], rjk)
        if mode == 2
            loggff += diploidfactor(s[j], h[j], p[j], pq[j], p̄[j])
        end
    end
    return exp(loggff)
end

function locuseffect(s, h, p, pq, m, p̄, r)
    q = 1 - p
    q̄ = 1 - p̄
    N = -s*h*(q̄ - q) + s*(1-2h)*(p̄*q - pq)
    D = m + r + s*(h - q + 2*(1-2h)*pq)
    # haploid -> sₑ = 2s, h=1/2 -> m + r + 2s(1/2 - (1-p)) = m + r - s(1 - 2p)
    return N/D
end

# when first generation of selection is a diploid migrant (SLiM?)
function diploidfactor(s, h, p, pq, p̄)
    q = 1 - p
    q̄ = 1 - p̄
    -s*(q̄ - q) + s*(1-2h)*(p̄*q̄ - pq)
end

function predict(s, h, u, Nₑ, p̄, mₑ)
    # p̄ = 0
    # p is locally beneficial, q is locally deleterious
    # migration introduces q
    q̄ = 1 - p̄
    d = Wright(-Nₑ*s, Nₑ*(mₑ*p̄ + u), Nₑ*(mₑ*q̄ + u), h)
    #@info -Nₑ*s, Nₑ*(mₑ*p̄ + u), Nₑ*(mₑ*q̄ + u), Nₑ, mean(d), mₑ
    # mean(d) is the frequency of the locally beneficial allele
    mean(d), expectedpq(d)
end

function gff(eqmodel::Equilibrium{M}, x) where M<:ZSF24
    @unpack model, Ep, Epq = eqmodel
    @unpack m, s, h, p̄, xs, rm = model
    loggff = 0.0
    for j=1:length(xs)
        rij = recrate(rm, x, xs[j])
        le = locuseffect(s[j], h[j], Ep[j], Epq[j], m, p̄[j], rij)
        loggff += le
        if mode == 2
            loggff += diploidfactor(s[j], h[j], Ep[j], Epq[j], p̄[j])
        end
    end
    exp(loggff)
end

function fixed_point_iteration(
        model::M, Ep, tol, nmax, α=1.0) where M<:ZSF24
    @unpack xs, Ne, m, s, h, u, p̄ = model
    it = 0
    Epq = Ep .* (1 .- Ep)
    while true
        gs = gff(model, Ep, Epq)
        ys = map(i->predict(s[i], h[i], u, Ne, p̄[i], m*gs[i]), 1:length(xs))
        _Ep = first.(ys)
        _Ep = (Ep .* (1-α)) .+ (α .* _Ep)
        if norm(Ep .- _Ep) < tol
            return Ep, Epq
        elseif it > nmax
            @warn "Not converged! # iterations > nmax (=$nmax)"
            @info "Retry with α = $(α/2)"
            fixed_point_iteration(model, Ep, tol, nmax, α/2)
        end
        _Epq = last.(ys)
        _Epq = (Epq .* (1-α)) .+ (α .* _Epq)
        Ep = _Ep
        Epq = _Epq
        it += 1
    end
end

