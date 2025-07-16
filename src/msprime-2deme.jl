
@with_kw struct TwoDemeCoalescent{T}
    L  :: Int
    r  :: T
    NA :: T  # source pop fwd in time
    NB :: T  # sink pop fwd in time
    m  :: T  # mig A->B fwd in time
    u  :: T = 0
    na :: Int
    nb :: Int
    ploidy :: Int = 2
end

function randts(θ::TwoDemeCoalescent; seed=rand(1:2^32))
    @unpack L, u, r, NA, NB, m, na, nb, ploidy = θ
    @assert m > 0.0
    msprime = pyimport("msprime")
    demo = msprime.Demography()
    demo.add_population(name="A", initial_size=NA)
    demo.add_population(name="B", initial_size=NB)
    demo.set_migration_rate(source="B", dest="A", rate=m)
    ts = msprime.sim_ancestry(
        samples=Dict("A"=> na, "B"=> nb), 
        recombination_rate=r, 
        sequence_length=L, 
        ploidy=ploidy, 
        random_seed=seed,
        demography=demo)
end

function randvars(θ::TwoDemeCoalescent, ts=randts(θ); seed=rand(1:2^32))
    @unpack u, na, nb, ploidy = θ
    msprime = pyimport("msprime")
    mts = msprime.sim_mutations(
        ts, rate=u, model=msprime.BinaryMutationModel(), 
        random_seed=seed
    )
    ns = (na + nb)*ploidy - 1
    xa = [i+1 for i=0:ns if mts.get_population(i) == 0]
    xb = [i+1 for i=0:ns if mts.get_population(i) == 1]
    idx = [xa ; xb]
    vars = map(mts.variants()) do var
        gs = var.genotypes.tolist()
        [ceil(Int, var.site.position); gs[idx]]
    end 
    colnames = ["POS" ; 
        ["A$i" for i=1:length(xa)] ; ["B$i" for i=1:length(xb)]]
    mts, DataFrame(permutedims(hcat(vars...)), colnames)
end

function coaltimes(mts, wins=[0, mts.sequence_length])
    sample_sets = [mts.samples(population=0), mts.samples(population=1)]
    ta = mts.diversity(mts.samples(population=0), mode="branch", windows=wins) / 2
    tb = mts.diversity(mts.samples(population=1), mode="branch", windows=wins) / 2
    tab = mts.divergence(sample_sets, mode="branch", windows=wins) / 2
    fst = mts.Fst(sample_sets, mode="branch", windows=wins)
    (ta=ta, tb=tb, tab=tab, fst=fst)
end

