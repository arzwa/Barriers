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
└   (E.Ep[1], E.Epq[1]) = (0.012914966757458318, 0.008115163867275742)

````

effective migration rates at map positions `0:0.1:1`

````julia
me_profile = map(y->Barriers.me(E, y), 0:0.1:1)
````

````
11-element Vector{Float64}:
 0.0023205422959092755
 0.0009565242006933221
 0.0012726071096870687
 0.000953147283045951
 0.0008468073320226461
 0.0011903483174785087
 0.0011763533781347687
 0.0007941370227800826
 0.0012601624010251766
 0.0024461195292587294
 0.0025948743151865132
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

