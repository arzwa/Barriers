

function mcmc(x, θ, proposals, n; every=10)
    @unpack ps, R, prior = θ
    L = length(ps)
    # row i contains the data to compute the gff at locus i
    # col j contains the contribution of locus j to all locus-specific gffs
    M = [i == j ? 0.0 : -x[j] * ps[j] / R[i,j] for i=1:L, j=1:L]
    ℓ = ℓhood(M, x, θ)
    π = sum(logpdf.(prior, x))
    x_ = copy(x)
    M_ = copy(M)
    X  = Matrix{Float64}(undef, n, L)
    ℓs = Vector{Float64}(undef, n)
    for k=1:n
        k % every == 0 && @info k, ℓ + π 
        for i=1:L
            x, M, ℓ, π = update!(i, x_, x, M_, M, ℓ, π, θ, proposals[i])
        end
        X[k,:] .= x
       ℓs[k]   = ℓ + π
    end
    return X, ℓs
end

function update!(i, x_, x, M_, M, ℓ, π, θ, prop)
    # update selection coefficient at locus i 
    x_[i], q = prop(x[i])
    M_[:,i] .= M[:,i] .* (x_[i]/x[i])
    # compute the likelihood, in principle this requires recomputing all
    ℓ_ = ℓhood(M_, x_, θ)
    π_ = π + logpdf(θ.prior, x_[i]) - logpdf(θ.prior, x[i])
    if log(rand()) <  ℓ_ + π_ - ℓ - π + q
        ℓ = ℓ_
        π = π_
        x[i] = x_[i]
        M[:,i] .= M_[:,i]
        accept!(prop)
    end
    return x, M, ℓ, π
end

# Do we need `Wright`? Or can we get away with the unnormalized pdf? I guess not:
# the normalization constant is only constant wrt to p, it is a function of the
# parameters...
function ℓhood(M, x, θ)
    @unpack ps, m, N, R, u = θ
    me = m*exp.(vec(sum(M, dims=2)))
    ds = Wright.(Ref(N), Ref(u), u .+ me, -x, Ref(0.0))
    ℓs = logpdf.(ds, ps)
    return sum(ℓs) 
end

function soften(ps)
    ϵ = minimum(filter(x->x!=0.0, ps))/10
    [p == 0.0 ? ϵ : p == 1 ? 1-ϵ : p for p in ps]   
end



# The same, but update a number of loci simultaneously
# ========================================================================
function mcmc1(x, θ, proposals, n; every=10)
    @unpack ps, R, prior, q = θ
    L = length(ps)
    # row i contains the data to compute the gff at locus i
    # col j contains the contribution of locus j to all locus-specific gffs
    M = [i == j ? 0.0 : -x[j] * ps[j] / R[i,j] for i=1:L, j=1:L]
    ℓ = ℓhood(M, x, θ)
    π = sum(logpdf.(prior, x))
    x_ = copy(x)
    M_ = copy(M)
    X  = Matrix{Float64}(undef, n, L)
    ℓs = Vector{Float64}(undef, n)
    for k=1:n
        k % every == 0 && @info k, ℓ + π 
        #K = min(L, rand(Geometric(q)+1))
        idxs = shuffle_partition(L, q)
        for idx in idxs
            x, M, ℓ, π = update1!(idx, x_, x, M_, M, ℓ, π, θ, proposals)
        end
        X[k,:] .= x
        ℓs[k]   = ℓ + π
    end
    return X, ℓs
end

shuffle_partition(N, k) = (collect ∘ partition)(shuffle(1:N), k)

function update1!(idx, x_, x, M_, M, ℓ, π, θ, proposals)
    # update selection coefficient at locus i 
    π_ = π 
    q_ = 0.0
    for i in idx
        x_[i], q = proposals[i](x[i])
        M_[:,i] .= M[:,i] .* (x_[i]/x[i])
        π_ += logpdf(θ.prior, x_[i]) - logpdf(θ.prior, x[i])
        q_ += q
    end
    # compute the likelihood, in principle this requires recomputing all
    ℓ_ = ℓhood(M_, x_, θ)
    if log(rand()) <  ℓ_ + π_ - ℓ - π + q_
        ℓ = ℓ_
        π = π_
        for i in idx
            x[i] = x_[i]
            M[:,i] .= M_[:,i]
            accept!(proposals[i])
        end
    end
    return x, M, ℓ, π
end


# Approximate approach for DFE inference
# =========================================================================
# Would be good to have a category of neutral loci: one more parameter p0

function egff(ps, R, α, θ)
    L  = length(ps)
    𝔼g = zeros(L)
    for i=1:L
        l𝔼gi = 0.0
        for j=1:L
            i == j && continue
            l𝔼gi += log(1 + θ*(ps[j]/R[i,j]))
        end
        𝔼g[i] = exp(-α * l𝔼gi)
    end
    return 𝔼g
end

function discretize(d, K)
    qstart = 1.0/2K
    qend = 1. - 1.0/2K
    xs = quantile.(d, qstart:(1/K):qend)
    xs .* (mean(d)*K/sum(xs))  # rescale by factor mean(d)/mean(xs)
end

function ℓhood2(g, x, θ)
    @unpack ps, m, N, R, u, K = θ
    L  = length(ps)
    me = m*g
    ℓ  = 0.0
    dfe = Gamma(x[1], x[2])
    ys  = discretize(dfe, K)
    for i=1:L
        ℓi = 0.0
        for k=1:K
            ℓi += pdf(Wright(N, u, u + me[i], -ys[k], 0.0), ps[i])
        end
        ℓ += log(ℓi) - log(K)
    end
    return ℓ
end

function mcmc2(x, θ, proposals, n; every=10)
    @unpack ps, R, prior = θ
    g = egff(ps, R, x[1], x[2])
    ℓ = ℓhood2(g, x, θ)
    π = logpdf(prior[1], x[1]) + logpdf(prior[2], x[2])
    x_ = copy(x)
    X  = Matrix{Float64}(undef, n, 2)
    ℓs = Vector{Float64}(undef, n)
    for k=1:n
        k % every == 0 && @printf "%6d %12.3f %9.3f %9.5f\n" k (ℓ + π) x[1] x[2] 
        for i=1:2
            x, ℓ, π = update2!(i, x_, x, ℓ, π, θ, proposals[i])
        end
        X[k,:] .= x
        ℓs[k]   = ℓ + π
    end
    return X, ℓs
end

function update2!(i, x_, x, ℓ, π, θ, prop)
    @unpack ps, R, prior = θ
    # update selection coefficient at locus i 
    x_[i], q = prop(x[i])
    # compute the likelihood, in principle this requires recomputing all
    g  = egff(ps, R, x_[1], x_[2])
    ℓ_ = ℓhood2(g, x_, θ)
    π_ = logpdf(prior[1], x_[1]) + logpdf(prior[2], x_[2])
    if log(rand()) <  ℓ_ + π_ - ℓ - π + q
        ℓ = ℓ_
        π = π_
        x[i] = x_[i]
        accept!(prop)
    end
    return x, ℓ, π
end

# Under this approach, it would be useful to classify loci according the DFE
# model, or better, sample from the posterior of their selection coefficient
# given the posterior for the DFE model.

# How do we do this? Consider we assume the posterior mean DFE. For a given
# locus i, we calculate the Egff. Then we know the distribution of allele
# frequencies for this site given some s (likelihood model), and we can use
# this to sample s. The problem is one-dimensional, so we can calculate the
# posterior by integration.


function posts(p, dfe, gff, θ; mins=1e-9, maxs=0.2)
    # P(s|p) = P(p|s)P(s)/P(p) = P(p|s)P(s)/∫P(p|s)P(s)ds
    # calculate the normalizing constant
    ℓ(s) = pdf(Wright(θ.N, θ.u, θ.u + θ.m*gff, -s, 0.0), p)
    Z, _ = quadgk(s->ℓ(s)*pdf(dfe,s), mins, maxs)
    post(s) = ℓ(s)*pdf(dfe,s) / Z
    return post 
end

function findquantile(f, q; xmin=1e-9, xmax=0.2)
    g(x) = quadgk(f, xmin, x)[1] - q
    find_zero(g, [xmin, xmax])
end

# with a class of neutral genes
# ===================================================================
# Could we resample which genes are neutral and which selected? It is something
# we would like to infer, and it need not necessarily make computation more
# tricky? Or does it? 
# P(selected=0|x) = P(x|selected=0)π0/(Σₛ P(x|selected=1,s)*P(s)*(1-π0) + P(x|selected=0)*π0)
function egff(ps, R, α, θ, π0)
    L  = length(ps)
    𝔼g = zeros(L)
    for i=1:L
        l𝔼gi = 0.0
        for j=1:L
            i == j && continue
            l𝔼gi += log(π0 + (1-π0)*(1 + θ*(ps[j]/R[i,j]))^(-α))
        end
        𝔼g[i] = exp(l𝔼gi)
    end
    return 𝔼g
end

function discretize(d, K)
    qstart = 1.0/2K
    qend = 1. - 1.0/2K
    xs = quantile.(d, qstart:(1/K):qend)
    xs .* (mean(d)*K/sum(xs))  # rescale by factor mean(d)/mean(xs)
end

function ℓhood3(g, x, θ; K=8)
    @unpack ps, m, N, R, u = θ
    L  = length(ps)
    me = m*g
    ℓ  = 0.0
    dfe = Gamma(x[1], x[2])
    π0  = x[3]
    ys  = discretize(dfe, K)
    for i=1:L
        ℓi = 0.0
        for k=1:K
            ℓi += pdf(Wright(N, u, u + me[i], -ys[k], 0.0), ps[i])
        end
        ℓi *= (1-π0)/K
        ℓi += π0*pdf(Wright(N, u, u + me[i], 0.0, 0.0), ps[i]) 
        ℓ  += log(ℓi)
    end
    return ℓ
end

function mcmc3(x, θ, proposals, n; every=10)
    @unpack ps, R, prior = θ
    g = egff(ps, R, x[1], x[2], x[3])
    ℓ = ℓhood3(g, x, θ)
    π = logpdf(prior[1], x[1]) + logpdf(prior[2], x[2]) + logpdf(prior[3], x[3])
    x_ = copy(x)
    X  = Matrix{Float64}(undef, n, 3)
    ℓs = Vector{Float64}(undef, n)
    for k=1:n
        k % every == 0 && @printf "%6d %12.3f %9.3f %9.5f %9.3f\n" k (ℓ + π) x[1] x[2] x[3] 
        for i=1:3
            x, ℓ, π = update3!(i, x_, x, ℓ, π, θ, proposals[i])
        end
        X[k,:] .= x
        ℓs[k]   = ℓ + π
    end
    return X, ℓs
end

function update3!(i, x_, x, ℓ, π, θ, prop)
    @unpack ps, R, prior = θ
    # update selection coefficient at locus i 
    x_[i], q = prop(x[i])
    # compute the likelihood, in principle this requires recomputing all
    g  = egff(ps, R, x_[1], x_[2], x_[3])
    ℓ_ = ℓhood3(g, x_, θ)
    π_ = logpdf(prior[1], x_[1]) + logpdf(prior[2], x_[2]) + logpdf(prior[3], x_[3])
    if log(rand()) <  ℓ_ + π_ - ℓ - π + q
        ℓ = ℓ_
        π = π_
        x[i] = x_[i]
        accept!(prop)
    end
    return x, ℓ, π
end
