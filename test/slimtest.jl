
# Single-locus mainland-island model
s = 0.04
m = 0.01
u = s/1000
N = 200
A  = Barriers.Architecture(
    [Barriers.DiploidLocus(2s, 0.5, u)], 
    zeros(1), zeros(1,1))
M  = Barriers.MainlandIslandModel(arch=A, m=m, N=N)
n  = 100_000

slim_script = Barriers.toslim(M, Barriers.linearmap(2,1.0), n, print_freqs=true)
path = "/home/arthur_z/dev/slim-sims/divsel-01.eidos"
write(path, slim_script)

buf = IOBuffer()
out = run(pipeline(`slim $path`, stdout=buf))
qss = split(String(take!(buf)), "\n")[end-n:end-1]

fun(x) = mapreduce(y->parse(Float64, y), +, split(x))
qs = map(fun, qss)

using WrightDistribution
d = Wright(2N*2s, 2N*(m + u), 2N*u, 0.5)
stephist(qs, norm=true, bins=0:0.01:1)
plot!(0:0.001:1, p->pdf(d, p))
