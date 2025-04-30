# 2024-11-25: modify sample3.jl to make use of the CL approach under the M
# model.
#
# The model we want is the following (subscript is means per window):
# - `m`  baseline migration
# - `ν`  density of selected sites
# - `s`  selection coefficient
# - `λ`  coalescence rate 
# - `α`  coalescence rate variation parameter (discrete Gamma)
# - `K`  K rate classes
# - `μᵢ` mutation rate for window i

# A more advanced model would involve an explicit model for BGS as well,
# instead of assuming ad hoc Nₑ variation... Note that we could fix μᵢ based on
# a global mutation rate (which sets the time scale for the other params)
# with/without additional local variation based on neutral substitution rates
# for instance.

struct CLSampler{T,V,X}
    priors    :: T  # NamedTuple with an entry for each parameter
    proposals :: V
    data      :: X
end

logprior(priors, θ) = mapreduce(x->logpdf(x...), +, zip(priors, θ))

struct CLData{V,T}
    X :: V          # site pattern counts for each window
    R :: Matrix{T}  # between window recrates
    L :: Vector{T}  # window lengths
    l :: Vector{T}  # relative window lengths
end
Base.length(d::CLData) = length(d.L)

function CLData(d::WindowedChromosome, X)
    n = length(d)
    R = window_recrates(d)
    L = [d[i][2][2] - d[i][2][1] for i=1:n]
    l = L ./ mean(L)
    L ./= sum(L)
    CLData(X, R, L, l)
end

function loglhood(data::CLData, G, θ)
    logmₑ = mepred(G, θ.m)
    loglhood(data.X, data.l, θ, logmₑ)
end

function loglhood(X::AbstractVector, l, θ, logmₑ)
    mapreduce(i->winlogpdf(X[i], l[i], θ, exp(logmₑ[i])), +, 1:length(X))
end

function winlogpdf(x, l, θ, m)
    @unpack K, λ, α, μ, γ = θ
    rs = discretize(Gamma(α, 1/α), K)
    ℓ = mapreduce(r->logpdfcl(m, μ, λ*r, λ*r, x, l*γ), +, rs) 
    ℓ - log(K)
end

function discretize(d, K)
    qstart = 1.0/2K
    qend = 1. - 1.0/2K
    xs = quantile.(d, qstart:(1/K):qend)
    xs .* (mean(d)*K/sum(xs))  # rescale by factor mean(d)/mean(xs)
end

function probs(G, θ, data::CLData, i, kmax=100, logeps=-12)
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
        p = logpdf(P, k-1) + loglhood(data.X, data.l, θ, logmₑ)
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

