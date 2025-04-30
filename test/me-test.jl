

# Compare against Sachdeva fig 1A.
L  = 80
s  = 0.02
Ns = 2.0
u  = s*0.005
ms = 0.0:0.02:0.8
Nₑ = ceil(Int64,Ns/s)

# haploid selection sₑ = 2s!
loci = fill(Barriers.DiploidLocus(2s, 0.5, u), L)
A  = Barriers.Architecture(loci, fill(Inf, L), fill(0.5, L, L))
M  = Barriers.MainlandIslandModel(arch=A, m=0.2s, N=Nₑ)
EM = Barriers.Equilibrium(M).Ep[1]

res = map(ms) do m
    M  = Barriers.MainlandIslandModel(arch=A, m=m*s, N=Nₑ)
    EM = Barriers.Equilibrium(M);
    EM.Ep[1]
end

plot(ms, res)

