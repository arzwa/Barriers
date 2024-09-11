# The goal is to devise a sampler.
# In the gIMble paper they estimates σ = s*ν and m.
# We could assume a fixed selection coefficient and infer the density ν of
# selected sites and background m.
# The model could be, for window i (windows iid)
#   m          ~ Exponential(m̄)
#   nᵢ         ~ Poisson(Lᵢ*ν)     # number of sel. sites in win i
#   xᵢ|nᵢ      ~ Dirichlet(nᵢ, α)  # stick breaking for sel. sites
# for the complete data
#   log(mₑ)|x  ~ Normal(log(mₑ[Aeschbacher model]), σ)  
# A sampler should iterate over windows?

"""
- `priors`: (s, ν, m, σ, α) 
"""
struct Sampler{T,V,X,W}
    priors    :: T  # NamedTuple with an entry for each parameter
    proposals :: V
    data      :: X
    meobs     :: W
end

struct State{U,T,V}
    θ :: U
    points  :: Vector{T}
    loci    :: Vector{T}  # locus positions in each window
    windowp :: Vector{V}  # posterior contribution of each window
    priorp  :: V          # global prior logpdf 
    lhood   :: V          # data loglikelihood
end

logprior(priors, θ) = mapreduce(x->logpdf(x...), +, zip(priors, θ))

function win_logpdf(θ, window, points)
    win, _ = window
    n  = length(points)
    l1 = logpdf(Poisson(θ.ν * length(win)), n) 
    n == 0 && return l1
    l2 = logpdf(Dirichlet(n, θ.α), points)
    return l1 + l2
end

function getwindowp(θ, data, points) 
    map(i->win_logpdf(θ, data[i], points[i]), 1:length(data))
end

function initialize(S::Sampler)
    @unpack priors, data, meobs = S
    θ = map(rand, priors)
    X = sample_windows(θ, data)
    points  = first.(X)
    loci    = getindex.(X, 2)
    windowp = last.(X)
    mepred = getme(θ, loci, data)
    priorp = logprior(priors, θ)
    lhood  = llikelihood(meobs, last.(mepred), θ.σ)
    State(θ, points, loci, windowp, priorp, lhood)
end

function sample_window(θ, data, i)
    win, (g1, g2) = data[i]
    d1 = Poisson(θ.ν * length(win))
    n  = rand(d1)
    if n == 0
        (Float64[], Float64[], logpdf(d1, n))
    else
        d2 = Dirichlet(n, θ.α)
        x  = rand(d2)
        z  = points2pos(x, g1, g2)
        (x, z, logpdf(d2, x) + logpdf(d1, n))
    end
end

sample_windows(θ, data) = map(i->sample_window(θ, data, i), 1:length(data))

function points2pos(xs, g1, g2)
    x = [0.0 ; cumsum(xs)]
    y = [(x[i] + x[i+1])/2 for i=1:length(x)-1]
    y .* (g2-g1) .+ g1
end

function getme(θ, loci, data)
    xs = vcat(loci...)
    model = AeschbacherModel(θ.m, fill(θ.s, length(xs)), xs, n=10)
    return me_profile(model, data)
end

function llikelihood(meobs, mepred, σ)
    sum(logpdf(Normal(0.0, σ), log.(meobs) .- log.(mepred)))
end

function mcmc(S::Sampler, X::State, n)
    @unpack priors, proposals, data, meobs = S
    for i=1:n
        # global parameters -- don't change m_e!
        for sym in keys(proposals)
            @unpack θ, priorp, windowp, loci, points, lhood = X
            proposal = getfield(proposals, sym)
            x_, q = proposal(getfield(θ, sym))
            θ_  = reconstruct(θ, sym=>x_)
            wp_ = getwindowp(θ_, data, points)  # XXX need to fix this
            pp_ = logprior(priors, θ_)
            ar  = pp_ + sum(wp_) - (priorp + sum(windowp)) + q
            if log(rand()) < ar
                X = State(θ_, points, loci, wp_, pp_, lhood)
                accept!(proposal)
            end
        end
        # iterate over windows
        #for j=1:length(data)
        #    @unpack θ, priorp, windowp, loci, lhood = X
        #    locij, lj = sample_window(θ, data, j)
        #    loci_ = copy(loci)
        #    loci_[j] = locij
        #    mepred = getme(θ, loci_, data)
        #    ll = llikelihood(meobs, last.(mepred), θ.σ)
        #    ar = ll + lj - (lhood + windowp[j])
        #    if log(rand()) < ar
        #        windowp[j] = lj
        #        X = State(θ, loci_, windowp, priorp, ll) 
        #    end
        #end
        @info X.θ
    end
end
