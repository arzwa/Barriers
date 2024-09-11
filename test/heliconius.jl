using CSV, DataFrames, Plots, Distributions, AdaptiveProposals
using Barriers

df = CSV.read("/home/arthur_z/dev/Barriers/data/64bpX500blocks-rho-cds-sumstats-me-fit-range-comb.tsv", DataFrame)
df18 = filter(x->x[:chrom] == 18, df)
plot(df18[:,:winMid], df18[:,"m_cyd..mel"], size=(700,200))
plot!(df18[1:5:end,:winMid], df18[1:5:end,"m_cyd..mel"])

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

L = 100
s = fill(-1e-4, L)
p = cumsum(rand(Dirichlet(L, 0.1))) .* Barriers.maplength(gmap)

model1 = AeschbacherModel(1e-6, s, p, n=20)
mes1 = me_profile(model1, d)
plot!(Barriers.plotcoords(mes1)...)

priors = (
    s=Dirac(-1e-4),
    ν=Exponential(50/chr_end),
    m=Exponential(1e-5),
    α=Dirac(1.0),
    σ=Dirac(1.0)
)

proposals = (
    ν=PositiveProposal(),
    m=PositiveProposal()
)

S = Barriers.Sampler(priors, proposals, d, data[:,"m_cyd..mel"])
X = Barriers.initialize(S)
@info length(vcat(X.loci...))
Barriers.mcmc(S, X, 1000)



