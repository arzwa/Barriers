module Calibration

export sims, calibrate_cl, CalibrationSims
using Barriers
using PyCall, StatsBase, Distributions, Parameters
using ProgressMeter, Combinatorics, Optim, QuadGK

"""
    CalibrationSims

A container for the parameters of the calibratoin simulations.
"""
@with_kw struct CalibrationSims{T}
    L  :: Int
    na :: Int
    nb :: Int
    m  :: T
    u  :: T
    r  :: T
    N  :: T
    α  :: T
    ploidy :: Int = 1
end

struct CalibrationSimResults{T,U,V}
    G   :: Matrix{Int}  # simulations
    Z   :: U            # pair data
    tbl :: Matrix{T}    # ABC table (distance, weight, m)
    prior :: V
end

function sims(θ::CalibrationSims;
        prior=Exponential(θ.m), nabc=10_000)
    @unpack na, nb, L = θ
    cmax = (na + nb)/(na*(na-1)*nb*(nb-1))
    G = msprime_sims(θ)
    A = 1:na; B=na+1:na+nb
    y = jsfs(G, A, B, L)
    Z, z = getpairdata(G, A, B, L)
    res = abc_m(y, prior, nabc, θ) 
    CalibrationSimResults(G, Z, res, prior)
end

function calibrate_cl(X::CalibrationSimResults, θ::CalibrationSims; 
        qabc=0.05, cmax=1.0)
    @unpack tbl, Z, prior = X
    z = sum(Z)
    q = quantile(tbl[:,1], qabc)
    pass = tbl[tbl[:,1] .< q,:]
    lnapprox = fit(LogNormal, pass[:,3])
    opt = Optim.optimize(c->obj(z, c, lnapprox, prior, θ), 0.0, cmax)
    copt = opt.minimizer
    pfun = clpost(z, copt, prior, θ)
    return copt, lnapprox, pfun
end

"""
    msprime_sims(θ::CalibrationSims)

Conduct simulations under the coalescent with recombination. Requires that
msprime is loaded.
```
msprime = pyimport("msprime")
```
"""
function msprime_sims(θ::CalibrationSims)
    @unpack L, u, r, N, α, m, na, nb, ploidy = θ
    @assert m > 0.0
    msprime = pyimport("msprime")
    demo = msprime.Demography()
    demo.add_population(name="A", initial_size=N*α)
    demo.add_population(name="B", initial_size=N)
    demo.set_migration_rate(source="A", dest="B", rate=m)
    ts = msprime.sim_ancestry(
        samples=Dict("A"=> na, "B"=> nb), 
        recombination_rate=r, 
        sequence_length=L, 
        ploidy=1, 
        demography=demo)
    mts = msprime.sim_mutations(
        ts, rate=u, model=msprime.BinaryMutationModel())
    ys = permutedims(
        hcat([var.genotypes.tolist() for var in mts.variants()]...))
end

"""
    jsfs(G, ia, ib, L, folded=true)

Obtain the joint site frequency spectrum, including invariant sites.
`G` is a genotype matrix, `ia` the indices of population A samples, `ib` for
population B samples. `L` is the sequence length (to count invariant sites).
"""
function jsfs(G::Matrix, ia, ib, L, folded=true)
    na = length(ia)
    nb = length(ib)
    sfs = zeros(na÷2+1,nb÷2+1)
    if size(G,2) == 0
        sfs[1,1] = 1.0 
        return sfs
    end
    ya = sum(G[:,ia], dims=2)
    yb = sum(G[:,ib], dims=2)
    if folded 
        ya = [y > na/2 ? na-y : y for y in ya]
        yb = [y > nb/2 ? nb-y : y for y in yb]
    end
    cs = countmap(collect(zip(ya, yb)))
    cs[(0,0)] = L - size(G,1)
    for i=0:na÷2, j=0:nb÷2
        haskey(cs, (i,j)) && (sfs[i+1,j+1] = cs[(i,j)]/L)
    end
    return sfs
end

#function jsfs_kldiv(ps, qs, eps=1e-10)
#    D = 0.0
#    k,l = size(ps)
#    for i=1:k, j=1:l
#        D += ps[i,j]*(log(ps[i,j] + eps) - log(qs[i,j] + eps))
#    end
#    return D
#end
#
#function jsfs_dist(y1, y2, eps=1e-10) 
#    jsfs_kldiv(y1, y2, eps) + jsfs_kldiv(y2, y1, eps)
#end

function getpairdata(sample, A, B, L)
    pairsA = combinations(A, 2)
    pairsB = combinations(B, 2)
    D = map(Iterators.product(pairsA, pairsB)) do (a, b)
        cols = [a ; b]
        Barriers.getcounts(sample[:,cols], L)
    end
    D, sum(D)
end

function abc_m(data, mproposal, nrep, θ)
    @unpack na, nb, L = θ
    @showprogress desc="Simulating $nrep replicates" map(1:nrep) do i
        m = rand(mproposal)
        w = pdf(mproposal, m)
        x = reconstruct(θ, m=m)
        G = msprime_sims(x)
        y = jsfs(G, 1:na, na+1:na+nb, L)
        d = sum((y .- data) .^ 2)
        [d, 1/w, m]
    end |> x->hcat(x...) |> permutedims
end

function clpost(z, c, mprior, θ)
    @unpack u, N, α = θ
    p(m) = Barriers.logpdfcl(m, u, 1/N, α/N, z, c) + logpdf(mprior, m)
    Z = quadgk(m->exp(p(m)), 0, Inf)[1]
    m->p(m) - log(Z)
end

# we optimize the KL divergence wrt `c`
function obj(z, c, ℓapprox, mprior, θ)
    pq = quantile(mprior, 0.999)
    p = clpost(z, c, mprior, θ)
    quadgk(m->exp(p(m))*(p(m) - logpdf(ℓapprox, m)), 0.0, pq)[1]
end

end
