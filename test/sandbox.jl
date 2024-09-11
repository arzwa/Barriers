using Barriers
using Barriers: WindowedChromosome, GeneticMap
using Barriers: MILocus, MIModel, AeschbacherModel, me_profile
using Distributions
using Plots

# a random genetic map
maplen = 1.0
α = 5.0
xs = [0.0 ; cumsum(rand(Dirichlet(100,α)) .* maplen)]
total = 20Mb
points = total ÷ (length(xs) - 1)
genetic_map = collect(zip(1:points:total, xs))
push!(genetic_map, (total, maplen))
gmap = GeneticMap(genetic_map)
d = WindowedChromosome(gmap, 200kb, 50kb) 

# L selected loci
L = 100
Ls = 0.7
xss = ceil.(Int64, cumsum(rand(Dirichlet(L+1, 10.))[1:L]) .* total)
pos = [d.geneticmap[x] for x in xss]
s = Ls/L
ss = -rand(Exponential(Ls/L), L)
m = s/2

# Aeschbacher model
model1 = AeschbacherModel(m, ss, pos)
mes1 = me_profile(model1, d)
ys1 = map(x -> ([extrema(x[1])...], [x[2], x[2]]), mes1)

# Diffusion model
Nes = 10.0
loci = [MILocus(ss[i], 0.5, s/100, 0.0, Nes/s) for i=1:L]
model2 = MIModel(m, pos, loci)
mes2 = me_profile(model2, d)
ys2 = map(x -> ([extrema(x[1])...], [x[2], x[2]]), mes2)

plot(first.(ys1), last.(ys1) ./ m, color=:black,
    legend=false, size=(700,200),
    ylabel="\$m_e\$", xlabel="physical position")
plot!(first.(ys2), last.(ys2) ./ m, color=:red, margin=3Plots.mm) 
xxs = 1:10kb:total
plot!(xxs, map(x->Barriers.me(model1, d.geneticmap[x]), xxs) ./ m, 
    color=:black, alpha=0.2)
plot!(xxs, map(x->Barriers.me(model2, d.geneticmap[x]), xxs) ./ m, 
    color=:red, alpha=0.2, margin=5Plots.mm)

