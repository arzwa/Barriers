# Retake on inference for the coarse-grained model
struct CLSampler{T,V,X}
    priors    :: T  # NamedTuple with an entry for each parameter
    proposals :: V
    data      :: X
end

@with_kw struct CLData{T}
    y :: Vector{Vector{Int}}      # counts for each window
    Δ :: Vector{T}                # window lengths in Morgan
    l :: Vector{T} = Δ ./ mean(Δ) # relative window length (for CL)
end
nwindow(y::CLData) = length(y.y)

struct State{T,V,M}
    θ     :: T 
    G     :: Matrix{V}
    logπ  :: V
    logℓ  :: V
    model :: M
end

function logstring(x::State, it, fields...)
    p = x.logπ + x.logℓ
    s = join([@sprintf("%.4e", getfield(x.θ, v)) for v in fields], " ")
    @sprintf "%7d %4e %s\n" it p s
end

function logprior(priors, θ) 
    mapreduce(x->logpdf(x...), +, zip(priors, θ))
end

function loglhood(data::CLData, G::Matrix, θ)
    logmₑ = mepred(G, θ.m)
    loglhood(data.y, data.l, θ, logmₑ)
end

function mepred(G, m) 
    vec(sum(G, dims=2)) .+ log(m)  # log mₑ
end

function loglhood(X::AbstractVector, l, θ, logmₑ)
    mapreduce(i->winlogpdf(X[i], l[i], θ, exp(logmₑ[i])), +, 1:length(X))
end

function winlogpdf(x, l, θ, m)
    @unpack α, λ, u, γ = θ
    logpdfcl(m, u, α*λ, λ, x, γ*l)
end

# XXX assuming complete divergence...
function log_gff_matrix(model::CoarseModel)
    n = length(model.X)
    G = zeros(n, n)
    [i == j ? log(gii(model, i)) : loggij(model, i, j) for i=1:n, j=1:n]
end

function initialize(S::CLSampler, θ=map(rand, S.priors))
    @unpack s, m, u = θ
    # [ν] is [# selected sites per map length]
    X = [rand(Poisson(θ.ν * l)) for l in S.data.Δ]
    model = CoarseModel(X=X, Δ=S.data.Δ, s=s, m=m, u=u, λ=0.0)  
    # XXX assuming complete divergence!
    G = log_gff_matrix(model)
    logπ = logprior(S.priors, θ) 
    logℓ = loglhood(S.data, G, θ)
    # this excludes P(Xᵢ|νLᵢ)...
    return State(θ, G, logπ, logℓ, model)
end

# XXX Neither ν_update nor X_update! involve MH steps -- we sample
# directly from the posterior
function ν_update(smplr, state)
    @unpack θ, G, model, logπ = state
    prior = smplr.priors.ν
    nsite = sum(model.X) / sum(model.Δ)
    posterior = Gamma(prior.α + nsite, 1/(1 + 1/prior.θ))
    ν_  = rand(posterior)
    @reset state.θ.ν = ν_ 
    @reset state.logπ = logprior(smplr.priors, state.θ) 
end

function X_update(smplr, state; windows=1:nwindow(smplr.data), kmax=100)
    state = deepcopy(state)
    for i in windows
        Xi_update!(state, smplr.data, i, kmax)
    end
    @reset state.logℓ = loglhood(smplr.data, state.G, state.θ)
end

function Xi_update!(state, data::CLData, i, kmax, logeps=-12)
    @unpack θ, G, model = state
    probs = Xi_posterior(state, data, i, kmax, logeps)
    # sample a new Xᵢ
    model.X[i] = sample(0:kmax, Weights(probs))
    # recalculate the gffs
    for j=1:size(G,1)
        G[j,i] = i == j ? log(gii(model, i)) : loggij(model, j, i)
    end
end

function Xi_posterior(state, data, i, kmax, logeps=-12)
    @unpack θ, G, model = state
    Po = Poisson(θ.ν*data.Δ[i])  # prior for number of selected sites in win
    g  = vec(sum(G, dims=2)) .- G[:,i]  
    gi = similar(g)
    # We will calculate the gffs & post. probs for Xᵢ ∈ [0,1,...kmax]
    # P(Xᵢ=k|data) = P(data|Xᵢ=k)P(Xᵢ=k)/∑ₖP(data|Xᵢ=k)P(Xᵢ=k)
    ps = fill(-Inf, kmax+1)
    k  = 0
    largest = -Inf
    while k <= kmax
        model.X[i] = k
        for j=1:length(g)
            gi[j] = i == j ? log(gii(model, i)) : loggij(model, j, i)
        end
        logmₑ = log(θ.m) .+ (g .+ gi)
        p = logpdf(Po, k) + loglhood(data.y, data.l, θ, logmₑ)
        ps[k+1] = p
        if p > largest 
            largest = p
        end
        k > 0 && (p < ps[k] && ps[k+1] < largest + logeps) && break
        k += 1
    end
    return lognormalize(ps)
end

function lognormalize(l)
   x = exp.(l .- maximum(l))
   return x ./ sum(x)
end

function mh_update(smplr, state, sym)
    @unpack logπ, logℓ = state
    @unpack priors, proposals, data = smplr
    proposal = getfield(proposals, sym)
    x_, q  = proposal(getfield(state.θ, sym))
    θ = reconstruct(state.θ, sym=>x_)  # not sure how to do this with Accessors
    model = hasfield(CoarseModel, sym) ? 
        reconstruct(state.model, sym=>x_) : state.model  
    # XXX assuming complete divergence!
    G = log_gff_matrix(model)
    logπ_ = logprior(priors, θ) 
    logℓ_ = loglhood(data, G, θ)
    ar  = logπ_ + logℓ_ - (logπ + logℓ) + q
    if log(rand()) < ar
        accept!(proposal)
        return State(θ, G, logπ_, logℓ_, model)
    else
        return state
    end
end

