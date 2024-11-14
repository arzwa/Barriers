using Pkg; Pkg.activate("/home/arthur_z/dev/Barriers")
using Barriers, StatsBase, Distributions, Parameters
using AdaptiveProposals, Plots, StatsPlots; plotsdefault();
include("test/heliconius-data.jl")

m = 1e-4
s = 2e-3
ν = 2e2  # density per M
σ = 10.
n = length(d)
priors = (
    m=Exponential(m), 
    s=Exponential(s), 
    ν=Dirac(ν))
#priors = (
#    m=Dirac(m), 
#    s=Dirac(s), 
#    ν=Dirac(ν))

# fake data
θ = (m=m, s=s, ν=ν, σ=σ)
data = Barriers.Data(d, zeros(n), σ)
X = [rand(Poisson(ν*data.L[i])) for i=1:n]
state = Barriers.initialize(θ, X, priors, data)
data  = Barriers.Data(d, exp.(log(m) .+ vec(sum(state.G, dims=2))), σ)
state = Barriers.initialize(θ, X, priors, data)

g = Barriers.predict_me(θ, X, data) ./ m
G = Barriers.g_matrix(θ, X, data)
gg = exp.(sum(G, dims=2))
all(g .≈ gg)

GG = copy(G)
Barriers.gcolumn!(@view(GG[:,1]), 1, 1, θ, data)

# makes sense? 
P = [Barriers.probs(state, data, i) for i=1:n]
X_ = map(P) do p
    sample(0:length(p)-1, Weights(p))
end
sum(X_), sum(X)

# mcmc
proposals = (
    ν=PositiveProposal(),
    s=PositiveProposal(),
    m=PositiveProposal()
)
S = Barriers.Sampler(priors, proposals, data)
res = Barriers.mcmc(S, state, 2000, m=20, kmax=10)

XX = hcat([x.X for x in res]...)
plot(mean(XX, dims=2))
plot!(X)

plot([x.θ.m for x in res])
hline!([m])

# Heliconius data
priors = (
    m=Gamma(0.1,2e-4), 
    s=Gamma(0.1,1e-3), 
    ν=Gamma(0.5,100.0))   # ν = 100 [/M] -> about 50 selected loci on chr 18
@info map(mean, priors)
@info map(x->quantile(x, [0.025, 0.975]), priors)

proposals = (
    ν=PositiveProposal(),
    s=PositiveProposal(),
    m=PositiveProposal()
)
θ     = map(mean, priors)
σ     = √(1/2)
data  = Barriers.Data(d, last.(meobs), σ)
X     = [rand(Poisson(θ.ν*data.L[i])) for i=1:length(data)]
state = Barriers.initialize(θ, X, priors, data)
S     = Barriers.Sampler(priors, proposals, data)

res_  = Barriers.mcmc(S, state, 110000, nw=50, kmax=100, every=10, thin=10)
res   = res_[1001:end]

zs = [last(d[i][1]) for i=1:length(d)]
XX = hcat([x.X for x in res]...)
ms = [x.θ.m for x in res]
ss = [x.θ.s for x in res]
vs = [x.θ.ν for x in res]
PP1 = plot(
    plot(log.(ms), color=:black, ylabel="\$\\log m\$"), 
    plot(log.(ss), color=:black, ylabel="\$\\log s\$"), 
    plot(log.(vs), color=:black, ylabel="\$\\log \\nu\$"), 
    size=(900,250), layout=(1,3), xlabel="iteration")

P1 = scatter(log.(ss), log.(vs), ms=1, color=:black, size=(500,500), 
    xlabel="\$\\log s\$", ylabel="\$\\log \\nu\$", legend=false)
P2 = scatter(log.(ss), log.(ms), ms=1, color=:black, size=(500,500), 
    xlabel="\$\\log s\$", ylabel="\$\\log m\$", legend=false)
P3 = scatter(log.(ms), log.(vs), ms=1, color=:black, size=(500,500), 
    xlabel="\$\\log m\$", ylabel="\$\\log \\nu\$", legend=false)
PP2 = plot(P1, P2, P3, size=(900,250), layout=(1,3), margin=4Plots.mm, legend=false)

plot(PP1, PP2, layout=(2,1), size=(900,500), legend=false)



bar(proportionmap(vec(sum(XX, dims=1))))
plot!(Poisson(mean(vs)))

function printtable(xs, d)
    println("| parameter | mean | 2.5% | 97.5% |")
    @printf "| --------- | ---- | ---: | ----: |\n"
    for (lab, xs) in zip(
            ["\$m\$", "\$s\$", "\$\\nu\$", "\$\\nu \\times s\\ [\\text{bp}^{-1}\$]"], 
            [ms, ss, vs, ss .* vs ./ Barriers.physlength(d)])
        xx = mean(xs)
        q1, q2 = quantile(xs, [0.025, 0.975])
        @printf "| %s | %.2e | %.2e | %.2e |\n" lab xx q1 q2
    end
end

printtable(xs, d)

# as expected these are anticorrelated.

p1=plot(priors.s, ylim=(0,1e3))
stephist!(ss, normalize=true)
p2=plot(priors.m, ylim=(0,1.2e6), xlim=(0,1e-5))
stephist!(ms, normalize=true)
p3=plot(priors.ν, ylim=(0,0.05))
stephist!(vs, normalize=true)
plot(p1,p2,p3,size=(800,250), layout=(1,3))

mm = hcat([exp.(Barriers.mepred(x.G, x.θ.m)) for x in res]...)
plot(zs, last.(meobs), color=:gray, ms=2, line=:steppre, 
    ylabel="\$m_e\$", xlabel="bp", bottom_margin=6Plots.mm, left_margin=3Plots.mm)
scatter!(zs, mean(mm, dims=2), size=(800,200), color=:firebrick, ms=1.5, 
    legend=false)
XX = hcat([x.X for x in res]...)
plot!(twinx(), zs, mean(XX, dims=2), color=:salmon, fill=true, alpha=0.5,
    line=:steppre, legend=false)

