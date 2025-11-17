# # Barriers
#
# This is actively developed research software. Everything will break from
# time to time.

using Barriers, Random
Random.seed!(112)

# Predict allele frequency divergence for $L$ unlinked loci
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

# With linkage, and variation in $s$
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

# effective migration rates at map position `y` (here middle of chromosome)
y = 0.5
Barriers.me(E, y)

# Model of Aeschbacher et al. 2017
AM = AeschbacherModel(m, [-A.loci[i].s*A.loci[i].h for i=1:L], x)
Barriers.me(AM, y)
