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
r = Barriers.recrate(genetic_map, 1, 2)

# Strongish selection
Ns = 50.0 
N  = 1000
s  = Ns/N
L  = 2
Ls = L*s

# Relatively weak migratiom
Nm = 1.0
m  = Nm/N
@info m/s

# Mutation is weak
u  = s/1000 

# Scattered across the genome, but not too close
Δ  = G÷L
xs = Δ÷2:Δ:G |> collect
ys = [genetic_map[x] for x in xs]

# Recombination rates between selected loci
R  = Barriers.recrates(genetic_map, xs)

# Mainland island model prediction
loci = fill(Barriers.DiploidLocus(s, 0.5, u), L)
A  = Barriers.Architecture(loci, ys, R)
M  = Barriers.MainlandIslandModel(arch=A, m=m, N=2N)
EM = Barriers.Equilibrium(M)

# Profile
ms = map(x->(x,Barriers.me(EM, x)), range(0, C, 1000))
fs = zip(first.(ms), 1 ./ (1 .+ 4N .* last.(ms))) |> collect
plot(fs, size=(900,200))
vline!(ys)

# Write a SLiM script
path = joinpath(@__DIR__, "slim", "test-01.eidos")

#slim_script = Barriers.toslim(
#    M, genetic_map, 10N, log=100, tsout="$path.ts", N1=N)
#write(path, slim_script)
out = run(pipeline(`slim $path`))

using PyCall
tskit = pyimport("tskit")
msprime = pyimport("msprime")

prefix = split(path, "-")[1]
pths = filter(x->startswith(x, prefix) && endswith(x, ".ts"), readdir(dirname(path), join=true))

function complete(ts, m, r)
    dem = msprime.Demography()
    dem.add_population(name="pop_0", initial_size=ts.num_individuals/2)
    dem.add_population(name="p1", initial_size=ts.num_individuals/2)
    dem.add_population(name="p2", initial_size=ts.num_individuals/2)
    dem.set_migration_rate(source="p2", dest="p1", rate=m)  # p1 -> p2 fwd in time
    tt = ts.dump_tables()
    tt.mutations.clear()
    ts = tt.tree_sequence()
    coalesced_ts = msprime.sim_ancestry(
        demography=dem, initial_state=ts, recombination_rate=r, ploidy=2)
end

function neutral_sims(ts, u, nrep)
    map(1:nrep) do _
        mts = msprime.sim_mutations(
            ts, rate=u, model=msprime.BinaryMutationModel())
    end
end

function fst_windows(ts, nwin)
    G = ts.sequence_length
    xw = map(x->ceil(Int64, x), range(0, G, nwin))
    smples = [ts.samples(population=p) for p in [1,2]]
    ts.Fst(smples, windows=xw, mode="branch")
end

res = map(pths) do pth
    @info pth
    ts = tskit.load(pth)
    ts = complete(ts, m, r)
    fst_windows(ts, 100)
end

loci = fill(Barriers.DiploidLocus(s, 0.5, u), L)
A  = Barriers.Architecture(loci, ys, R)
M  = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
EM = Barriers.Equilibrium(M)
ms = map(x->(x,Barriers.me(EM, x)), range(0, C, 1000))
fs = zip(first.(ms), 1 ./ (1 .+ 8N .* last.(ms))) |> collect
plot(fs, size=(900,200))
vline!(ys)
xw = map(x->genetic_map[ceil(Int64, x)], range(1, G, 100))[2:end]
plot!(xw, mean(res), line=:steppre)

model = Barriers.AeschbacherModel(m, -s, ys)
ams = map(x->m .* Barriers.gff(model, x), 0:0.001:1)
afs = 1 ./ (1 .+ 8N .* ams)
plot!(0:0.001:1, afs)
