
@with_kw struct BPModel{T,V,R<:RecombinationMap,N}
    m  :: T
    s  :: Vector{T}
    xs :: Vector{V}
    rm :: R
    Ne :: N
    u  :: T
end

struct Equilibrium{M,T}
    model :: M
    Ep    :: T
    Epq   :: T
end

function Equilibrium(model; 
        p0=ones(length(model.xs)),
        tol=1e-9, nmax=100, α=1.)
    Ep, Epq = fixed_point_iteration(model, p0, tol, nmax, α)
    Equilibrium(model, Ep, Epq)
end

# Selective gff
gff(eq::Equilibrium) = gff(eq.model, eq.Ep)
function gff(model::BPModel, Ep::AbstractVector) 
    @unpack s, rm, xs = model
    gff(rm, xs, s .* Ep)
end

function gff(rm::RecombinationMap, xs::AbstractVector, sp::AbstractVector)
    L = length(xs)
    map(1:L) do k
        gff(rm, xs[k], xs, sp, sp[k], k)
    end
end

# Neutral gff
function gff(eq::Equilibrium{M}, x) where M<:BPModel
    @unpack s, rm, xs = eq.model
    gff(rm, x, xs, eq.Ep .* s)
end

me(eq::Equilibrium, x) = eq.model.m*gff(eq, x)

function gff(rm::RecombinationMap, x, xs, sp, spk=0., k=0)
    loggff = 0.0
    j = 1
    if k == 0 
        k = findlast(z->z < x, xs)
        k = isnothing(k) ? 0 : k
        j = 0
    end
    S = spk   # 0 for neutral gff
    for i=k-j:-1:1
        r = recrate(rm, xs[i], x)
        a = sp[i]/(r + S)
        loggff -= log1p(a)
        S += sp[i]
    end
    S = spk
    for i=k+1:length(xs)
        r = recrate(rm, xs[i], x)
        a = sp[i]/(r + S)
        loggff -= log1p(a)
        S += sp[i]
    end
    return exp(loggff)
end

function fixed_point_iteration(
        model::BPModel, Ep, tol, nmax, α=1.0)
    @unpack xs, Ne, m, s, u = model
    it = 0
    while true
        gs = gff(model, Ep)
        _Ep = map(1:length(xs)) do i
            d = Wright(-2Ne*s[i], Ne*u, Ne*(m*gs[i] + u), 0.5)
            mean(d)
        end
        _Ep = (Ep .* (1-α)) .+ (α .* _Ep)
        if norm(Ep .- _Ep) < tol
            return _Ep, fill(NaN, length(_Ep))
        elseif it > nmax
            @warn "Not converged! # iterations > nmax (=$nmax)"
            @info "Retry with α = $(α/2), nmax = $(2nmax)"
            fixed_point_iteration(model, Ep, tol, 2nmax, α/2)
        end
        Ep = _Ep
        it += 1
    end
end

"""
    bc_fitnesses

Calculate expected backcross fitnesses and frequencies to 'zeroth order',
i.e. an F1 carries on average half of the migrants deleterious alleles, a
BC1 half of that, etc.
"""
function bc_fitnesses(eq, kmax)
    @unpack m, s, xs = eq.model
    ws = map(0:kmax-1) do k
        exp(-sum([s[i]*eq.Ep[i] for i=1:length(xs)])/2^k)
    end 
    fs = m*cumprod(2 .* ws)
    fs, ws
end

"""
    bc_ancestries

Calculate expected proportion of migrant ancestry in backcross generations.
That is, hs[k] is the expected proportion of migrant ancestry in the whole
population that is in BCk individuals. ∑ₖhs[k] is an approximation for the
total proportion of migrant ancestry in the population at equilibrium.
"""
function bc_ancestries(eq, kmax)
    @unpack m, s, xs, rm = eq.model
    R = Fwd.rec_matrix(rm, xs)
    _, ws = bc_fitnesses(eq, kmax)
    S = sum(eq.Ep .* s)
    hs = zeros(kmax)
    hs[1] = m*ws[1]
    for k=1:kmax-1
        hs[k+1] = exp(-S*avr(R, k))*hs[k]
    end
    return hs
end

function avr(R, t)
    RT = (1 .- R) .^ t
    L = size(R,1)
    r̄ = 0.
    for i=2:L
        for j=1:i-1
            r̄ += RT[i,j]
        end
    end
    r̄ *= 2/(L*(L-1))
end
