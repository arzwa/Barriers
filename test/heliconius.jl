using CSV, DataFrames, Plots, Distributions, AdaptiveProposals, Parameters
using Barriers
plotsdefault();

df = CSV.read("/home/arthur_z/dev/Barriers/data/64bpX500blocks-rho-cds-sumstats-me-fit-range-comb.tsv", DataFrame)
df18 = filter(x->x[:chrom] == 18, df)
plot(df18[:,:winMid], df18[:,"m_cyd..mel"], size=(700,200))
plot!(df18[1:5:end,:winMid], df18[1:5:end,"m_cyd..mel"])

map(3:7) do i
    plot(df18[1:5:end,:winMid], df18[1:5:end,i], title="$(names(df18)[i])")
end |> x->plot(x..., layout=(5,1), size=(300,700), margin=3Plots.mm)

plot(df18[:,:winMid], df18[:,3])
plot!(df18[:,:winMid], df18[:,4])
plot!(df18[:,:winMid], df18[:,5], ls=:dot, size=(700,200))

# this seems to be how they calculated fst
# i.e. (πB - πW)/(πB + πW)
fst = (df18[:,5] .- (df18[:,3] .+ df18[:,4]) ./ 2) ./ (df18[:,5] .+ (df18[:,3] .+ df18[:,:4]) ./2)

lmap = CSV.read("/home/arthur_z/dev/Barriers/data/linkage_map.cydno.cm.tsv", DataFrame)
lmap = filter(x->x[:Chromosome] == "chr18", lmap)
lmap = lmap[.!nonunique(lmap[:,[:cM]]),:]
gmap = GeneticMap(collect(zip(lmap[:,:End] .+ 1, lmap[:,:cM] ./ 100)))

# take this as the end of the chromosome (a window is approximately 100kb) 
chr_end = maximum(df18[:,:winMid])+50kb

# issue!
chr_end > Barriers.physlength(gmap)

# we should extrapolate based on the average recombination rate?
gmap = Barriers.extrapolate_to(gmap, chr_end)
chr_end == Barriers.physlength(gmap)

# nonoverlapping windows
data = df18[1:5:end,:]
xs = data[:,:winMid]
zs = [1 ; floor.(Int64,[(xs[i] + xs[i-1])/2 for i=2:length(xs)]) ; chr_end+1]
wins = [zs[i]:zs[i+1]-1 for i=1:length(zs)-1]

# the object
d = WindowedChromosome(gmap, wins)

meobs = [win => m for (win, m) in zip(first.(d.data), data[:,"m_cyd..mel"])]
plot(Barriers.plotcoords(meobs)..., size=(700,200))


# -----------------------
L = 100
s = fill(-1e-4, L)
p = cumsum(rand(Dirichlet(L, 0.1))) .* Barriers.maplength(gmap)

model1 = AeschbacherModel(1e-6, s, p, n=20)
mes1 = me_profile(model1, d)
plot!(Barriers.plotcoords(mes1)...)


# should first find the globally best fitting ν assuming equally spaced loci,
# use this to set the prior and initialize.
priors = (
    s=Exponential(1e-4),
#    s=Dirac(1e-4),
    ν=Exponential(100/chr_end),
    m=Exponential(1e-6),
    α=Dirac(1.0),
    σ=Dirac(0.15)
)
proposals = (
    ν=PositiveProposal(),
    s=PositiveProposal(),
    m=PositiveProposal()
)
S = Barriers.Sampler(priors, proposals, d, data[:,"m_cyd..mel"])
X = Barriers.initialize(S)
res = Barriers.mcmc(S, X, 2500, m=50, l=3)


res = [res ; Barriers.mcmc(S, res[end], 3000, m=10, l=3)]

res_ = res[501:end]
Y = permutedims(hcat(map(x->[x.θ.m, x.θ.s, x.θ.ν], res_)...))
plot(plot(Y[:,1]), plot(Y[:,2]), plot(Y[:,3]),
    size=(1000,200), layout=(1,3))
map(zip(eachcol(Y), ["m", "s", "ν"])) do (y, v)
    @printf("%s, %4.2e, (%4.2e, %4.2e)\n", 
        v, mean(y), quantile(y, [0.025, 0.975])...)
end;

xx = map(x-> Barriers.plotcoords(x.mepred), res_)
avg = vec(mean(hcat(last.(xx)...), dims=2))
winsize = [length(x[1]) for x in d.data]
nd = map(x->length.(x.loci) ./ winsize, res_)
nd = mean(nd)
n = map(x->length.(x.loci), res_)

xs, ys = Barriers.plotcoords(meobs)
P1 = plot(xs, ys, color=:gray, label="gIMble", legend=:topright,
    ylabel="\$m_e\$", 
    title="Aeschbacher et al. model")
plot!(xs, avg, label="fit", color=:firebrick)
P2 = plot(xs, repeat(nd, inner=2), label="", color=:lightgray, 
    ylabel="local \$\\nu\$", fill=true, xlabel="location (Mb)")
plot(P1, P2, layout=(2,1), size=(800,350), guidefont=9, 
    xlim=(0,chr_end), xticks=(0:5Mb:chr_end, (0:5Mb:chr_end) ./ Mb),
     margin=4Plots.mm)


vv = map(x->x.θ.ν * x.θ.s, res_)
mm = map(x->x.θ.m, res_)
@printf("νs = %4.2e (%4.2e, %4.2e)\n", mean(vv), quantile(vv, [0.025,0.975])...)
@printf("m  = %4.2e (%4.2e, %4.2e)\n", mean(mm), quantile(mm, [0.025,0.975])...)



pp1 = histogram(Y[:,2], norm=true)
plot!(priors.s)
pp2 = histogram(Y[:,1], norm=true)
plot!(priors.m)
pp3 = histogram(Y[:,3], norm=true)
plot!(priors.ν)
plot(pp1,pp2,pp3,layout=(1,3), size=(700,200))
# I don't understand how we can have so much information about these params
# independently?

m, s, ν = mean(Y, dims=1)
priors = (
    s=Dirac(s),
    ν=Dirac(ν),
    m=Dirac(m),
    α=Dirac(1.0),
    σ=Dirac(0.5)
)
proposals = ()
S = Barriers.Sampler(priors, proposals, d, data[:,"m_cyd..mel"])
X = Barriers.initialize(S)
resb = Barriers.mcmc(S, X, 2500, m=50, l=3)

resb = [resb; Barriers.mcmc(S, resb[end], 50, m=50, l=3)]

xx = map(x-> Barriers.plotcoords(x.mepred), resb)
avg = vec(mean(hcat(last.(xx)...), dims=2))
winsize = [length(x[1]) for x in d.data]
nd = map(x->length.(x.loci) ./ winsize, resb)
nd = mean(nd)
n = map(x->length.(x.loci), resb)

xs, ys = Barriers.plotcoords(meobs)
P1 = plot(xs, ys, color=:gray, label="gIMble", legend=:topright,
    ylabel="\$m_e\$", 
    title="Aeschbacher et al. model")
plot!(xs, avg, label="fit", color=:firebrick)
P2 = plot(xs, repeat(nd, inner=2), label="", color=:lightgray, 
    ylabel="local \$\\nu\$", fill=true, xlabel="location (Mb)")
plot(P1, P2, layout=(2,1), size=(800,350), guidefont=9, 
    xlim=(0,chr_end), xticks=(0:5Mb:chr_end, (0:5Mb:chr_end) ./ Mb),
     margin=4Plots.mm)


# local density
# -------------------------------------------------------------------------

priors = (
    s=Exponential(1e-4),
    νσ=Dirac(0.5),
    ν=Normal(log(100/chr_end), 0.5),
    m=Exponential(1e-6),
    α=Dirac(1.0),
    σ=Dirac(0.3)
)
proposals = (
    ν=AdaptiveProposal(kernel=Normal(), trans=identity),
    s=PositiveProposal(),
    m=PositiveProposal(),
    νs=[AdaptiveProposal(kernel=Normal(), trans=identity) for i=1:length(d)]
)

S = Barriers.Sampler(priors, proposals, d, data[:,"m_cyd..mel"])
X = Barriers.initialize2(S)

@info length(vcat(X.loci...))
res = Barriers.mcmc2(S, X, 500, m=5, l=3)


# distributed
# -------------------------------------------------------------------------
using Distributed
using Plots, StatsPlots; plotsdefault()
addprocs(4)
@everywhere using Pkg; 
@everywhere Pkg.activate("/home/arthur_z/dev/Barriers/");
@everywhere using Barriers, AdaptiveProposals, Distributions, DataFrames, CSV

ress = pmap(1:4) do _
    priors = (
        s=Exponential(1e-4),
    #    s=Dirac(1e-4),
        ν=Exponential(100/chr_end),
        m=Exponential(1e-6),
        α=Dirac(1.0),
        σ=Dirac(0.10)
    )
    proposals = (
        ν=PositiveProposal(),
        s=PositiveProposal(),
        m=PositiveProposal()
    )
    S = Barriers.Sampler(priors, proposals, d, data[:,"m_cyd..mel"])
    X = Barriers.initialize(S)
    res = Barriers.mcmc(S, X, 1000, m=50, l=4)
    res = [res ; Barriers.mcmc(S, X, 4000, m=10, l=4)]
end

Ys = map(ress) do res
    res_ = res[1:end]
    Y = permutedims(hcat(map(x->[x.θ.m, x.θ.s, x.θ.ν], res_)...))
end

map(1:3) do j
    plot(); [plot!(Ys[i][:,j]) for i=1:3]; plot!()
end |> x->plot(x..., layout=(1,3), size=(900,300))
# Having issues converging...

res_ = ress[1][1001:end]
xx = map(res_) do x
    Barriers.plotcoords(x.mepred)
end
avg = vec(mean(hcat(last.(xx)...), dims=2))
winsize = [length(x[1]) for x in d.data]
nd = map(x->length.(x.loci) ./ winsize, res_)
nd = mean(nd)

xs, ys = Barriers.plotcoords(meobs)
P1 = plot(xs, ys, color=:gray, label="gIMble", legend=:topright,
    ylabel="\$m_e\$")
plot!(xs, avg, label="fit")
P2 = plot(xs, repeat(nd, inner=2), label="", color=:lightgray, 
    ylabel="local \$\\nu\$", fill=true, xlabel="location (Mb)")
plot(P1, P2, layout=(2,1), size=(800,350), guidefont=9, 
    xlim=(0,chr_end), xticks=(0:5Mb:chr_end, (0:5Mb:chr_end) ./ Mb),
     margin=4Plots.mm)
