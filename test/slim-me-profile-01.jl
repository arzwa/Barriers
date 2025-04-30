# Compare different mₑ approximations
using Pkg; Pkg.activate("/home/arthur_z/dev/Barriers")
using Barriers
using DataFrames, CSV, StatsBase, Distributions
using Plots; plotsdefault()
using Random

# A 10Mb chromosome, covering one Morgan
G = 10Mb
C = 2.0 
genetic_map = Barriers.linearmap(G, C)
@info Barriers.recrate(genetic_map, 1, G)
r = Barriers.recrate(genetic_map, 1, 2)

# Strongish selection
Ns = 20.0 
N  = 1000
s  = Ns/N
L  = 40
Ls = L*s

# Relatively weak migratiom
Nm = 0.1
m  = Nm/N
@info m/s

# Mutation is weak
u  = s/200 

# Scattered across the genome, but not too close
rng = Random.seed!(123)
α  = 1.0
xs = cumsum(rand(rng, Dirichlet(L, α)))
xs = [((i == 1 ? 0.0 : xs[i-1]) + xs[i])/2 for i=1:L]
xs = map(x->ceil(Int64, x), xs .* G)
ys = [genetic_map[x] for x in xs]

# Recombination rates between selected loci
R  = Barriers.recrates(genetic_map, xs)

# Mainland island model prediction
loci = fill(Barriers.DiploidLocus(2s, 0.5, u), L)
A  = Barriers.Architecture(loci, ys, R)
M  = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
EM = Barriers.Equilibrium(M)

# Profile
ms = map(x->(x,Barriers.me(EM, x)), range(0, C, 1000))
plot(ms, size=(900,200))
vline!(ys)

# Write a SLiM script
path = "/home/arthur_z/dev/slim-sims/divsel-03.eidos"
slim_script = Barriers.toslim(M, genetic_map, 10N, log=100, tsout="$path.ts", N1=N)
write(path, slim_script)
out = run(pipeline(`slim $path`))

using PyCall
tskit = pyimport("tskit")
mts = tskit.load("/home/arthur_z/dev/Barriers/data/ts.2024-12-13.trees")

muts = Barriers.ts_get_mutations(mts)
data = Barriers.get_counts_windows(muts, 100kb) 
plot(permutedims(hcat(data./ 100kb...)[2:end,:]))

# fit without architecture model, once this works we can worry about fitting
# the architecture
using Turing

@model function tmodel3(data, μ=5e-8)
    n = length(data)
    τ ~ Exponential(1)
    m ~ Normal(5.5, 5.0)
    l ~ Normal(5.5, 5.0)
    lms ~ MvNormal(fill(m, n), τ)
    ms = exp.(lms)
    for i=1:length(data)
        ps = Barriers.probs(ms[i]*μ, μ, exp(l)*μ, exp(l)*μ)
        if !isprobvec(ps)
            Turing.@addlogprob! -Inf
        else
            data[i] ~ Multinomial(sum(data[i]), ps) 
        end
    end
end

chn = sample(tmodel3(data), NUTS(), 500)

MM = exp.(chn.value[:,4:length(data)+4-1]) .* 2e-9
LL = exp.(chn.value[:,length(data)+4:2length(data)+4-1]) .* 2e-9

function plotinf(x, X, qq1=0.025, qq2=0.975; color2=:firebrick, kwargs...)
    yy = vec(mean(X, dims=1))
    q1 = map(x->quantile(x, qq1), eachcol(X))
    q2 = map(x->quantile(x, qq2), eachcol(X))
    scatter(x, yy, color=color2, ms=2)
    map(1:length(yy)) do i
        plot!([x[i], x[i]], [q1[i], q2[i]], label="", 
            color=color2, lw=1.5, alpha=0.5)
    end
    plot!(legend=false)
end

plotinf(50kb:100kb:G, 1 ./ LL)

plotinf(50kb:100kb:G, MM)






using AdaptiveProposals
winL = fill(Barriers.distance(genetic_map, 1, 100kb), length(data)) 
winl = winL ./ mean(winL)
R    = Barriers.recrates(genetic_map, collect(100kb÷2:100kb:G))

cldata = Barriers.CLData(data, R, winL, winl)

priors = (
    K=Dirac(1), 
    μ=Dirac(5e-8), 
    λ=Exponential(1/2000), 
    α=Dirac(1.0),  
    γ=Dirac(0.05),
    #m=Gamma(1,1e-4), 
    m=Dirac(0.0001),
    #s=Gamma(1,0.02), 
    s=Dirac(s),
    ν=Gamma(1,40.0)
)

proposals = (
    s=PositiveProposal(),
    m=PositiveProposal(),
    λ=PositiveProposal(),
#    α=PositiveProposal()
)

θ     = map(mean, priors)
Y     = [rand(Poisson(θ.ν*cldata.L[i])) for i=1:length(cldata)]
state = Barriers.initialize(θ, Y, priors, cldata)
S     = Barriers.CLSampler(priors, proposals, cldata)

res  = Barriers.mcmc(S, state, 5000, nw=10, kmax=20, every=5, thin=1)

res  = Barriers.mcmc(S, res[end], 20000, nw=100, kmax=20, every=100, thin=2)

res=res[1000:end]
map([:m, :ν, :s, :λ]) do sym
    xs = map(x->getfield(x.θ, sym), res)
    if sym==:λ
        xs = 1 ./ 2xs
    end
    plot(xs, color=:black, ylabel="\$\\log $sym\$")
end |> x->plot(x..., size=(900,400))

zs = 50kb:100kb:G-50kb-1
mm = hcat([exp.(Barriers.mepred(x.G, x.θ.m)) for x in res]...)
qq = map(x->quantile(x, [0.025, 0.975]), eachrow(mm)) 
q1 = first.(qq)
q2 = last.(qq)
Em = mean(mm, dims=2)
P1 = plot(zs, Em, ribbon=(Em .- q1, q2 .- Em),  lines=:steppre, label="95% interval")
plot!([genetic_map[z] for z in first.(ms)], last.(ms), color=:black, line=:steppre, 
    ylabel="\$m_e\$", xlabel="bp", size=(800,200), 
    label="true \$m_e\$", legend=:topleft,
    bottom_margin=6Plots.mm, left_margin=3Plots.mm, right_margin=3Plots.mm)

plot!(twinx(), zs, mean(XX, dims=2), color=:salmon, fill=true, alpha=0.5,
    legend=:topright, label="№ selected sites",
    line=:steppre)


