
struct Sampler{T,V,X,W}
    priors    :: T  # NamedTuple with an entry for each parameter
    proposals :: V
    data      :: X
    meobs     :: W
end

# ν for each window
struct LocalDensityState{T,V,W}
    θ      :: T           # global parameters
    νs     :: Vector{V}   # local density of selected sites
    ℓs     :: Vector{V}   # log prior density contribution per window
    points :: Vector{W}   # Dirichlet points in each window
    loci   :: Vector{W}   # loci locations in each window
    mepred :: Vector{V}   # predicted mₑ profile
    priorp :: V           # global log prior density
    lhood  :: V           # log likelihood
end

function mcmc2(S::Sampler, X::LocalDensityState, n; m=30, l=5)
    @unpack priors, proposals, data, meobs = S
    map(1:n) do k
        @printf("%6d, %5.3e, %5.3e, %5.3e, %9.5f\n", 
            k, X.θ.m, X.θ.ν, X.θ.s, X.lhood + X.priorp)
        # global parameters -- ν doesn't change m_e
        # m changes m_e in a very simple way.
        for sym in keys(proposals)
            sym == :νs && continue  # HACKy
            @unpack θ, νs, ℓs, points, loci, mepred, priorp, lhood = X
            proposal = getfield(proposals, sym)
            x_, q = proposal(getfield(θ, sym))
            θ_  = reconstruct(θ, sym=>x_)
            ℓs_ = window_logpdf(θ_, νs, points, data)
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
            ar  = ll_ + pp_ + sum(ℓs_) - 
               (lhood + priorp + sum(ℓs)) + q  # stupid to do sum(windowp)...
            if log(rand()) < ar
                X = LocalDensityState(θ_, νs, ℓs_, points, loci, mepred_, pp_, ll_)
                accept!(proposal)
            end
        end
        # iterate over windows
        for j=rand(1:length(data), m)
            X = window_update!(X, S, j, l)
        end
        X
    end
end

function window_update!(X, S, j, l)
    @unpack priors, proposals, data, meobs = S
    @unpack θ, νs, ℓs, points, loci, priorp, lhood, mepred = X
    old = (loci[j], νs[j])
    prop = proposals.νs[j]
    νj, q = prop(νs[j])
    pointsj, locij = sample_window(νj, θ.α, data[j])
    ℓj = window_logpdf(θ, νj, pointsj, length(data[j][1])) 
    loci[j] = locij
    mepred_ = getme!(θ, loci, data, copy(mepred), j, n=l)
    ll = llikelihood(meobs, mepred_, θ.σ)
    ar = ll + ℓj - (lhood + ℓs[j]) + q
    if log(rand()) < ar
        νs[j] = νj
        ℓs[j] = ℓj
        points[j] = pointsj
        X = LocalDensityState(θ, νs, ℓs, points, loci, mepred_, priorp, ll) 
        accept!(prop)
    else
        loci[j] = old[1]
    end
    return X
end

function initialize2(S::Sampler)
    @unpack priors, data, meobs = S
    θ       = map(rand, priors)
    νs      = rand(Normal(θ.ν, θ.νσ), length(data))
    X       = sample_windows(θ, νs, data)
    points  = first.(X)
    loci    = last.(X)
    windowp = window_logpdf(θ, νs, points, data)
    mepred  = getme(θ, loci, data)
    priorp  = logprior(priors, θ)
    lhood   = llikelihood(meobs, mepred, θ.σ)
    LocalDensityState(θ, νs, windowp, points, loci, mepred, priorp, lhood)
end

function sample_window(ν, α, window)
    win, (g1, g2) = window
    n = rand(Poisson(exp(ν) * length(win)))
    n == 0 && return (Float64[], Float64[])
    x = rand(Dirichlet(n, α))
    z = points2pos(x, g1, g2)
    return (x, z)
end

sample_windows(θ, νs, data) = map(
    i->sample_window(νs[i], θ.α, data[i]), 1:length(data))

function window_logpdf(θ, νs::Vector, points, data) 
    map(i->window_logpdf(θ, νs[i], points[i], 
        length(data[i][1])), 1:length(data))
end

function window_logpdf(θ, ν::T, points, winlen) where T<:Real
    n  = length(points)
    l1 = logpdf(Poisson(exp(ν) * winlen), n) + logpdf(Normal(θ.ν, θ.νσ), ν)
    n == 0 && return l1
    l2 = logpdf(Dirichlet(n, θ.α), points)
    return l1 + l2
end

function logprior(priors, θ) 
    mapreduce(x->logpdf(x...), +, zip(priors, θ))
end

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
