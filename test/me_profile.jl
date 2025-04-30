# Compare different mₑ approximations
using Pkg; Pkg.activate("/home/arzwa/dev/Barriers")
using Barriers
using DataFrames, CSV, StatsBase, Distributions
using Plots; plotsdefault()
using Random

# A 10Mb chromosome, covering one Morgan
G = 10Mb
C = 1.0
genetic_map = Barriers.linearmap(G, C)
@info Barriers.recrate(genetic_map, 1, G)
r = Barriers.recrate(genetic_map, 1, 2)

# Strongish selection
Ns = 20.0 
N  = 1000
s  = Ns/N
L  = 10
Ls = L*s

# Relatively weak migratiom
Nm = 1.0
m  = Nm/N
@info m/s

# Mutation is weak
u  = s/1000 

# Scattered across the genome, but not too close
Δ  = G ÷ L
xs = (Δ÷2):Δ:G |> collect
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
#plot(ms, size=(900,200))
plot(first.(ms), 1 ./ (1 .+ last.(ms) .* 4N))
vline!(ys)

# Write a SLiM script
path = "/home/arthur_z/dev/slim-sims/divsel-04.eidos"
slim_script = Barriers.toslim(
    M, genetic_map, 10N, log=100, tsout="$path.ts", N1=N)
write(path, slim_script)
out = run(pipeline(`slim $path`))

using PyCall
tskit = pyimport("tskit")
msprime = pyimport("msprime")
ts = tskit.load("$path.ts").simplify()
N = 1000
A = ts.samples()[1:2*N]
B = ts.samples()[2*N+1:end]

# First we have to complete the forward simulation backwards in time
dem = msprime.Demography()
dem.add_population(name="p1", initial_size=ts.num_individuals/2)
dem.add_population(name="p2", initial_size=ts.num_individuals/2)
dem.set_migration_rate(source="p2", dest="p1", rate=m)  # p1 -> p2 fwd in time

# remove mutations (obscure error due to mutation times, guess it has to do
# with the initialization of the mainland)
nt = ts.dump_tables()
nt.mutations.clear()
nt = nt.tree_sequence()

coalesced_ts = msprime.sim_ancestry(
        demography=dem,
        initial_state=nt,
        recombination_rate=r,
        ploidy=2)

function sample_sets(ts, each=2)
    [sample(ts.samples(population=p.id), each, replace=false) for p in ts.populations()]
end

μ = 1e-7
L = ceil(Int64, coalesced_ts.sequence_length)
w = map(x->round(Int64, x), range(0, L-1, length=101))
w[end] += 1  # something went wrong
res = map(1:50) do _
    mts = msprime.sim_mutations(
        coalesced_ts, rate=μ, model=msprime.BinaryMutationModel())
    f = mts.Fst(sample_sets(mts), windows=w, mode="branch")
end

plot([genetic_map[x+1] for x in w[1:end-1]], mean(res))
plot!(first.(ms), 1 ./ (1 .+ last.(ms) .* 4N))


