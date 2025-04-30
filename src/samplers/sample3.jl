# 2024-11-13: This is a nice sampler which uses a rather crude mₑ approximation
# and exploits conjugacy to sample ν.
struct Sampler{T,V,X}
    priors    :: T  # NamedTuple with an entry for each parameter
    proposals :: V
    data      :: X
end

logprior(priors, θ) = mapreduce(x->logpdf(x...), +, zip(priors, θ))

struct Data{V,T}
    ℓ :: V          # sampling distribution (MvNormal)
    R :: Matrix{T}  # between window recrates
    L :: Vector{T}  # window lengths
end
Base.length(d::Data) = length(d.L)

function Data(d::WindowedChromosome, meobs, σ)
    n = length(d)
    D = product_distribution([Normal(log(meobs[i]), σ) for i=1:n])
    R = window_recrates(d)
    L = [d[i][2][2] - d[i][2][1] for i=1:n]
    L ./= sum(L)
    Data(D, R, L)
end

function Data(d::Vector, σ)
    nelem = length(d)
    xs = map(d) do (dd, meobs)
        n = length(dd)
        D = [Normal(log(meobs[i]), σ) for i=1:n]
        R = window_recrates(dd)
        L = [dd[i][2][2] - dd[i][2][1] for i=1:n]
        D, R, L
    end
    D = product_distribution(vcat(getindex.(xs, 1)...))
    L = vcat(getindex.(xs, 3)...)
    L ./= sum(L)
    T = length(L)  # total number of windows
    R = fill(0.5, T, T)
    x0 = 1
    for Ri in getindex.(xs, 2)
        a, _ = size(Ri)
        R[x0:x0+a-1,x0:x0+a-1] .= Ri
        x0 += a 
    end
    Data(D, R, L)
end

struct State{T,U}
    θ      :: T
    X      :: Vector{Int64}
    G      :: Matrix{U}
    latent :: Vector{U}
    priorp :: U 
    lhood  :: U
end

function initialize(θ, X, priors, data)
    G      = g_matrix(θ, X, data)
    lhood  = loglhood(data, G, θ)
    latent = latent_logpdf(θ.ν, X, data)
    priorp = logprior(priors, θ)
    State(θ, X, G, latent, priorp, lhood)
end

"""
    g_matrix

This computes the effect of each window on each other window on the log gff,
and stores it in a matrix. `G[i,i]` is the contribution of selected loci within
window `i` on the average gff in window `i`. `G[i,j]` is the contribution of
window `j` to the average gff in window `i`. We have `sum(G[i,:]) == log(gᵢ)`.

Importantly, when the number of selected loci in window `j` changes, this only
affects `G[:,j]`, allowing for efficient updates.
"""
g_matrix(θ, X, data) = g_matrix(θ, X, data.L, data.R)

function g_matrix(θ, X, L, R)
    @unpack s, m = θ
    n = length(L)
    G = zeros(n,n)
    for i=1:n, j=1:n
        if i == j
            G[i,i] = local_lg(X[i], L[i], s, m)  # log me
        else
            G[i,j] = -X[j]*s/(m + s + R[i,j]) 
        end
    end
    return G
end

# Assumes equally spaced loci in the window, and integrate over the window to
# get average g.
function local_lg(X, L, s, m)
    X == 0 && return 0.0
    d  = L/X
    x0 = d/2
    xs = x0:d:L
    ḡ, _ = quadgk(x->exp(_lg(xs, s, m, x)), 0.0, L)
    # XXX opportunity for speed-up? FastGaussQuadrature?
    return log(ḡ / L)
end

function _lg(xs, s, m, xi)
    lgi = 0.0
    for xj in xs
        rij = recrate(abs(xj - xi)) 
        lgi -= s/(m + s + rij)
    end
    return lgi
end

mepred(G, m) = vec(sum(G, dims=2)) .+ log(m)  # log mₑ

function latent_logpdf(ν, X, data)
    [logpdf(Poisson(ν*data.L[i]), X[i]) for i=1:length(X)]
end

function loglhood(data::Data, G, θ)
    logmₑ = mepred(G, θ.m)
    logpdf(data.ℓ, logmₑ)
end

# assuming complete divergence
#function predict_me(θ, X, data::Data)
#    [predict_me(θ, X, data, i) for i=1:length(data)]
#end
#
#function predict_me(θ, X, data::Data, i)
#    @unpack s, m = θ
#    @unpack R, L = data
#    lgi = local_lg(X[i], L[i], s, m)  # log me
#    for j=1:length(data)
#        (j == i || X[j] == 0) && continue
#        lgi -= X[j]*s/(m + s + R[i,j]) 
#    end
#    return m*exp(lgi)
#end

function lognormalize(l)
   x = exp.(l .- maximum(l))
   return x ./ sum(x)
end

function probs(G, θ, data, i, kmax=100, logeps=-12)
    @unpack ν, m, s = θ
    L = data.L[i] 
    P = Poisson(ν*L)
    g = vec(sum(G, dims=2)) .- G[:,i]
    gi = similar(g)
    ps = Float64[]
    k  = 1
    largest = -Inf
    while k < kmax
        gcolumn!(gi, i, k-1, data.L, data.R, s, m)
        logmₑ = log(θ.m) .+ (g .+ gi)
        p = logpdf(P, k-1) + logpdf(data.ℓ, logmₑ)
        push!(ps, p)
        if p > largest 
            largest = p
        end
        if k > 1
            (p < ps[k-1] && ps[k] < largest + logeps) &&  break
        end
        k += 1
    end
    lognormalize(ps)
end

function gcolumn!(g, i, k, L, R, s, m) 
    for j=1:length(g)
        g[j] = (j == i) ? local_lg(k, L[i], s, m) : -k*s/(m + s + R[i,j]) 
    end    
end

function gcolumn!(g, i, k, θ, data)
    @unpack L, R = data
    gcolumn!(g, i, k, L, R, θ.s, θ.m)
end

function ν_update(state, data, prior::Gamma, priors)
    @unpack θ, G, X, lhood = state
    posterior = Gamma(prior.α + sum(X), 1/(1 + 1/prior.θ))
    ν_  = rand(posterior)
    θ_  = @set θ.ν=ν_
    pp_ = logprior(priors, θ_)
    latent_ = latent_logpdf(θ_.ν, X, data)
    state = State(θ_, X, G, latent_, pp_, lhood)
end

function mh_update(state, data, priors, proposal, sym)
    @unpack θ, X, G, priorp, lhood, latent = state
    x_, q  = proposal(getfield(θ, sym))
    θ_  = reconstruct(θ, sym=>x_)  # not sure how to do this with Accessors
    pp_ = logprior(priors, θ_)
    G_  = g_matrix(θ_, X, data)
    #latent_ = copy(latent)
    ll_ = loglhood(data, G_, θ_)
    ar  = ll_ + pp_ - (lhood + priorp) + q
    @info "update" θ_.m θ.m ll_ lhood pp_ priorp q
    if log(rand()) < ar
        accept!(proposal)
        return State(θ_, X, G_, latent, pp_, ll_)
    else
        return state
    end
end

function mcmc(S, state::V, n; 
        nw=30, every=10, kmax=20, thin=10) where V<:State
    @unpack priors, proposals, data = S
    res = Vector{V}(undef, n÷thin)
    try
        for it=1:n
            if (it-1) % every == 0 
                printlog(state, proposals, it-1)
            end
            # global parameters
            if typeof(priors.ν) <: Gamma
                state = ν_update(state, data, priors.ν, priors)
            end
            for sym in keys(proposals)
                state = mh_update(state, data, priors, getfield(proposals, sym), sym)
            end
            # iterate over windows
            X = copy(state.X)
            G = copy(state.G)
            l = copy(state.latent)
            for j=rand(1:length(data), nw)
                p = probs(G, state.θ, data, j, kmax)
                k = sample(0:length(p)-1, Weights(p))
                X[j] = k
                l[j] = logpdf(Poisson(state.θ.ν * data.L[j]), k)  
                # XXX what do we store `l` for? not needed...
                gcolumn!(@view(G[:,j]), j, k, state.θ, data)
            end
            lhood  = loglhood(data, G, state.θ)
            state  = State(state.θ, X, G, l, state.priorp, lhood)
            if it % thin == 0 
                res[it÷thin] = state
            end
        end
    catch e
        @error "Interrupted" e
        return res[filter(i->isdefined(res,i), 1:length(res))]
    end
    return res
end

function printlog(state, proposals, it)
    @printf "%6d " it
    for f in [keys(proposals)...; :ν]
        x = getfield(state.θ, f)
        @printf "%s=%5.3e " f x
    end
    print("\n")
end


