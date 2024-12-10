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

struct State{U,T,V,W}
    θ       :: U
    points  :: Vector{T}
    loci    :: Vector{T}  # locus positions in each window
    windowp :: Vector{V}  # posterior contribution of each window
    priorp  :: V          # global prior logpdf 
    lhood   :: V          # data loglikelihood
    mepred  :: Vector{W}
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

function initialize(S::Sampler, θ=map(rand, S.priors))
    @unpack priors, data, meobs = S
    X = sample_windows(θ, data)
    points  = first.(X)
    loci    = getindex.(X, 2)
    windowp = last.(X)
    mepred  = getme(θ, loci, data)
    priorp  = logprior(priors, θ)
    lhood   = llikelihood(meobs, mepred, θ.σ)
    State(θ, points, loci, windowp, priorp, lhood, mepred)
end

function sample_window(ν, α, window)
    win, (g1, g2) = window
    d1 = Poisson(ν * length(win))
    n  = rand(d1)
    if n == 0
        (Float64[], Float64[], logpdf(d1, n))
    else
        d2 = Dirichlet(n, α)
        x  = rand(d2)
        z  = points2pos(x, g1, g2)
        (x, z, logpdf(d2, x) + logpdf(d1, n))
    end
end

sample_windows(θ, data) = map(
    i->sample_window(θ.ν, θ.α, data[i]), 1:length(data))

function points2pos(xs, g1, g2)
    x = [0.0 ; cumsum(xs)]
    y = [(x[i] + x[i+1])/2 for i=1:length(x)-1]
    y .* (g2-g1) .+ g1
end

# XXX can be computationally much improved if we only recompute locally under
# the finite-MSP model...
function getme(θ, loci, data; n=2)
    xs = vcat(loci...)
    model = AeschbacherModel(θ.m, -θ.s, xs, n=n)
    return me_profile_vec(model, data)
end

function getme!(θ, loci, data, mes, j; n=2)
    xs = vcat(loci...)
    model = AeschbacherModel(θ.m, -θ.s, xs, n=n)
    return me_profile_vec_local!(model, data, mes, loci, j)
end

function llikelihood(meobs, mepred, σ)
    sum(logpdf(Normal(0.0, σ), log.(meobs) .- log.(mepred)))
end

"""
- `n`: number of iterations
- `m`: number of windows to resample per iteration
- `l`: number of MSPs to consider on both sides of a selected locus
"""
function mcmc(S::Sampler, X::State, n; m=30, l=5, every=10)
    @unpack priors, proposals, data, meobs = S
    map(1:n) do k
        k % every == 0 && (@printf("%6d, n%4d, m%5.3e, ν%5.3e, s%5.3e, π%9.5f\n", 
            k, length(vcat(X.loci...)), X.θ.m, X.θ.ν, X.θ.s, 
            X.lhood + X.priorp + sum(X.windowp)))
        # global parameters -- ν doesn't change m_e
        # m changes m_e in a very simple way.
        for sym in keys(proposals)
            @unpack θ, priorp, windowp, loci, points, lhood, mepred = X
            proposal = getfield(proposals, sym)
            x_, q = proposal(getfield(θ, sym))
            θ_  = reconstruct(θ, sym=>x_)
            wp_ = getwindowp(θ_, data, points) 
            pp_ = logprior(priors, θ_)
            if sym == :m
                mepred_ = (θ_.m/θ.m) .* mepred
                ll_ = llikelihood(meobs, mepred_, θ.σ)
            elseif sym == :s
                mepred_ = getme(θ_, loci, data, n=l)
                ll_ = llikelihood(meobs, mepred_, θ.σ)
            else
                mepred_ = copy(mepred)
                ll_ = lhood
            end
            ar  = ll_ + pp_ + sum(wp_) - 
               (lhood + priorp + sum(windowp)) + q  # stupid to do sum(windowp)...
            if log(rand()) < ar
                X = State(θ_, points, loci, wp_, pp_, ll_, mepred_)
                accept!(proposal)
            end
        end
        # iterate over windows
        for j=rand(1:length(data), m)
            X = window_update!(X.θ.ν, X, S, j, l)
        end
        X
    end
end

function window_update!(ν::Real, X, S, j, l)
    @unpack priors, proposals, data, meobs = S
    @unpack θ, priorp, windowp, points, loci, lhood, mepred = X
    old = loci[j]
    pointsij, locij, lj = sample_window(ν, θ.α, data[j])
    loci[j] = locij
    mepred_ = getme!(θ, loci, data, copy(mepred), j, n=l)
    ll = llikelihood(meobs, mepred_, θ.σ)
    ar = ll + lj - (lhood + windowp[j])
    if log(rand()) < ar
        windowp = copy(windowp)
        points = copy(points)
        loci = copy(loci)
        windowp[j] = lj
        points[j] = pointsij
        X = State(θ, points, loci, windowp, priorp, ll, mepred_) 
    else
        loci[j] = old
    end
    return X
end



