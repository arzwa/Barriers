using Sewall, AdaptiveProposals, Distributions

function soften(ps)
   ϵ = minimum(filter(x->x!=0.0, ps))/10
   [p == 0.0 ? ϵ : p == 1 ? 1-ϵ : p for p in ps]
end

function egff(ps, R, α, θ)
    L  = length(ps)
    𝔼g = zeros(L)
    for i=1:L
        l𝔼gi = 0.0
        for j=1:L
            i == j && continue
            l𝔼gi += log(1 + θ*(ps[j]/R[i,j]))
            # this works only for weak recombination, because the mgf hack
            # for dealing with Gamma distributed s does not work if s also
            # appears in the denominator in the exponent. i.e. we are using the
            # gamma mgf to compute E[exp(-s*p/r)] but this does not work for
            # E[exp(-s*p/(r + m + s*p)]
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
#            ℓi += pdf(Wright(N, u, u + me[i], -ys[k], 0.0), ps[i])
            ℓi += pdf(Wright(-N*ys[k], N*u, N*(u + me[i]), 0.5), ps[i])
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

function posts(p, dfe, gff, θ; mins=1e-9, maxs=0.2)
    # P(s|p) = P(p|s)P(s)/P(p) = P(p|s)P(s)/∫P(p|s)P(s)ds
    # calculate the normalizing constant
    @info (θ.N, θ.u, θ.u + θ.m*gff, -maxs, 0.0)
    @unpack N, u, m = θ
    ℓ(s) = pdf(Wright(-N*s, N*u, N*(u + m*gff), 0.5), p)
    Z, _ = quadgk(s->ℓ(s)*pdf(dfe,s), mins, maxs)
    post(s) = ℓ(s)*pdf(dfe,s) / Z
    return post
end

function findquantile(f, q; xmin=1e-9, xmax=0.2)
    g(x) = quadgk(f, xmin, x)[1] - q
    find_zero(g, [xmin, xmax])
end

# Simulate a large thing
rng = Random.seed!(321)
maplength = 4.0
Nes= 50
T  = 250
L  = 100
Ls = 0.7
s̄  = Ls/L
Ne = Nes/s̄
u  = s̄/200
k  = 1
N  = Ne2N(Ne, k)
ps = L/T
α  = 10.0
xs  = cumsum(rand(rng, Dirichlet(fill(α, T)))) .* maplength
xs .-= xs[1]
ds  = [xs[i] - xs[i-1] for i=2:T]  # pairwise distances
rs  = Sewall.haldane.(ds)
selidx = sort(sample(rng, 1:T, L, replace=false))
neuidx = setdiff(1:T, selidx)
loci = [Locus(0.0, 0.0, 0.0, u, u) for i=1:T]
κ   = 4/1
dfe = Gamma(κ, s̄/κ)
sss = rand(rng, dfe, L)
sss .*= (s̄/mean(sss))
for (i,j) in enumerate(selidx)
    loci[j] = Locus(-sss[i], 0.0, 0.0, u, u)
end
A = Architecture(loci, rs)
# mainland allele frequencies at mut-sel balance, neutral loci at 0.5
y = [x.s1 == 0.0 ? 0.5 : max(min(0.9999, 1-(-u/x.s1)), 0.0001) for x in loci]
mainland = HWLEMainland(y)
D = Deme(N=N, k=k, A=A)
m = 0.5*s̄
M = MainlandIsland(D, m, 0.0, mainland)
p = mean.(equilibrium(M))
q = 1 .- p;
x_, y_ = Sewall.meprofile(M, 500)
g = Sewall.eqgff(M)

P1 = plot(x_, y_, size=(900,200))
plot!(twinx(), x_, 1 ./ (1 .+ 2Ne .* y_), size=(900,200), legend=false, color=2)

rng = Random.seed!(323)
P = initpop(rng, M, zeros(T))
ff(P, _, _) = vec(mean(P.diploids, dims=1))
pop, res = simulate!(rng, M, P, 20000, ff)

X = hcat(res...)
ps = soften(X[:,end])

P2 = sticks(xs[neuidx], ps[neuidx], color=:black, legend=false)
sticks!(xs[selidx], ps[selidx], color=:firebrick, legend=false)
sticks!(xs[selidx], - q[selidx], ylim=(-1,1))
sticks!(xs[neuidx], - q[neuidx], color=:gray, ylim=(-1,1))
plot(P1, P2, layout=(2,1), size=(1200,400))


#θ = (N=M.deme.Ne, m=M.mhap, R=A.R, ps=1 .- ps, u=A[1].u01,
#    prior=Exponential(0.02), q=10)
#proposals = [PositiveProposal(kernel=Normal(0,0.1)) for i=1:T]
#Y, l = mcmc1(rand(θ.prior, T) ./ 10, θ, proposals, 10000, every=2)
prior = [Exponential(1), Exponential(0.02)]
proposals = [PositiveProposal(kernel=Normal(0,0.1)) for i=1:2]
θ = (N=M.deme.Ne, m=M.mhap, R=A.R, ps=1 .- ps, u=A[1].u01, prior=prior, K=10)
Y1, l1 = mcmc2(rand.(θ.prior), θ, proposals, 1500, every=10)


L = length(ps)
a, b = mean(Y1, dims=1)
gs  = egff(θ.ps, θ.R, a, b)
dfe = Gamma(a,b)
maxs = 100.0 / θ.N
pps = map(1:L) do i
    fp = posts(θ.ps[i], dfe, gs[i], θ, maxs=maxs)
    Es, _ = quadgk(s->s*fp(s), 0, maxs)
    q1 = findquantile(fp, 0.05, xmax=maxs)
    q2 = findquantile(fp, 0.95, xmax=maxs)
    @info i, Es, q1, q2
    Es, q1, q2, fp
end

Es = getindex.(pps, 1)
truth = [-A[i].s1 for i=1:length(A)]
P0 = sticks(xs, -Es, color=:firebrick, size=(1200,200),
    legend=false)
sticks!(xs, truth, color=:black)
plot!(x_, y_ ./ m, color=:black, legend=false)
plot!(xs, gs)

q1 = getindex.(pps, 2)
q2 = getindex.(pps, 3)
selected = selidx
# observed data
P1 = sticks(xs, θ.ps, legend=false, ms=3, alpha=0.2, color=:black,
    xlabel="", framestyle=:default,
    title="\$L = $(length(selected))\$\n",
    titlefont=12,
    ylabel="",
    ylim=(0,1))

# inference compared to truth
P2 = plot(framestyle=:default, xshowaxis=false, legend=:bottomright)
smin = maximum(truth)
sticks!(xs, width=3, truth, size=(1200,300), color=4, ylim=(-smin,smin),
    label="true \$s\$", ylabel="\$s\$")
sticks!(xs, -Es, width=3, color=2, label="\$\\mathbb{E}[s|p]\$")
hline!([0], color=:lightgray, label="")
scatter!(xs, zeros(L), color=:black, ms=1, label="", xticks=false)
# Posterior uncertainty
yts = (-4:1:-0.5, map(x->"\$10^{$(ceil(Int, x))}\$", -4:1:-0.5))
q1 .= max.(Ref(1.3e-4), q1)
P3 = plot(framestyle=:default, ylim=(-4., -0.5), xlabel="", legend=false, yticks=yts)
map(i->plot!(fill(xs[i],2), log10.([q1[i], q2[i]]), color=:lightgray, lw=2,
    alpha=1, label=""), setdiff(1:L, selected))
map(i->plot!(fill(xs[i],2), log10.([q1[i], q2[i]]), color=:firebrick,
    lw=2, alpha=1, label=""), selected)
scatter!(xs, log10.(truth), color=:black, size=(900,200), ms=2)
#    annotate!(5, -11.5, text("\$\\mathrm{Map\\ position}\$", :center, 8))
plot(P1, P2, P3, P0, layout=grid(4,1,heights=[0.1,0.3,0.3,0.3]), size=(900,600))

