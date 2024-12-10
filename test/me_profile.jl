# Compare different mₑ approximations
using Pkg; Pkg.activate("/home/arthur_z/dev/Barriers")
using Barriers
using DataFrames, CSV, StatsBase, Distributions
using Plots; plotsdefault()
using Random

# A 10Mb chromosome, covering one Morgan
G = 10Mb
C = 1.0 
genetic_map = Barriers.linearmap(G, C)
@info Barriers.recrate(genetic_map, 1, G)

# Strongish selection
Ns = 50.0 
N  = 1000
s  = Ns/N
Ls = 1.0
L  = ceil(Int64, Ls/s)

# Relatively weak migratiom
Nm = 0.5
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

loci = fill(Barriers.DiploidLocus(2s, 0.5, u), L)
A  = Barriers.Architecture(loci, ys, R)
M  = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
EM = Barriers.Equilibrium(M)

ms = map(x->(x,Barriers.me(EM, x)), range(0, C, 1000))
plot(ms, size=(900,200))
vline!(ys)




