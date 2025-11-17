```@meta
EditURL = "doc/README.jl"
```

# Barriers

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

With linkage, and variation in $s$

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
└   (E.Ep[1], E.Epq[1]) = (0.9375400966168393, 0.05576030226430959)

````

effective migration rates at map position `y`

````julia
y = 0.12
Barriers.me(E, y)
````

````
0.0006400073198132228
````

Model of Aeschbacher et al. 2017

````julia
AM = AeschbacherModel(m, [-A.loci[i].s*A.loci[i].h for i=1:L], x)
Barriers.me(AM, y)
````

````
6.453165767309359e-5
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

