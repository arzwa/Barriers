```@meta
EditURL = "doc/README.jl"
```

Barriers

This is actively developed research software. Everything will break from
time to time.

````julia
using Barriers, Random
````

Predict allele frequency divergence for $L$ unlinked loci

````julia
L = 50
s = 0.02    # selection coeff.
h = 0.5     # dominance coeff.
m = s/2     # migration rate
u = s*1e-3  # mutation rate
N = 500     # haploid population size  (no. haploid genomes)
l = [Barriers.DiploidLocus(2s, 0.5, u) for i=1:L]
R = fill(0.5,L,L)  # recombination rate matrix
x = fill(NaN, L)   # genomic locations, use NaN for unlinked architecture
A = Barriers.Architecture(l, fill(NaN, L), fill(0.5,L,L))
M = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
E = Barriers.Equilibrium(M)
@info "Beneficial allele freq., heterozygosity " E.Ep[1], E.Epq[1]
````

````
┌ Info: Beneficial allele freq., heterozygosity 
└   (E.Ep[1], E.Epq[1]) = (0.8995028953305044, 0.08496170472012593)

````

With linkage, and variation in s

````julia
L = 25
s = 0.02     # selection coeff.
h = 0.5      # dominance coeff.
m = s/2      # migration rate
u = s*1e-3   # mutation rate
N = 500      # haploid population size  (no. haploid genomes)
l = [Barriers.DiploidLocus(2s*randexp(), 0.5, u) for i=1:L]
x = sort(rand(L))  # genomic locations (1M chromosome)
A = Barriers.Architecture(l, x)
M = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
E = Barriers.Equilibrium(M)
@info "Beneficial allele freq., heterozygosity " E.Ep[1], E.Epq[1]
````

````
┌ Info: Beneficial allele freq., heterozygosity 
└   (E.Ep[1], E.Epq[1]) = (0.8656379736355391, 0.10932751995806925)

````

effective migration rates at map positions 0:0.1:1

````julia
me_profile = map(y->Barriers.me(E, y), 0:0.1:1)
````

````
11-element Vector{Float64}:
 0.0021560892425549335
 0.0010355591792402085
 0.0013441986913604984
 0.0014542946482371927
 0.0013913870808799208
 0.0006520226819396586
 0.00033754825712916905
 0.000925726179664812
 0.0015075011446559533
 0.00158032633833215
 0.0014547634902671402
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

