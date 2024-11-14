using Barriers
using Test

using Barriers: MILocus

m = 0.001
s = 0.005
u = s/1000
N = 200/s
locus = Barriers.MILocus(2s, 0.5, u, 0.0, N) 
Barriers.predict(locus, m)

mss = 0:0.01:2
plot(mss, map(x->Barriers.predict(locus, s*x)[1], mss))
vline!([1])

Ls = 1.0
L  = ceil(Int, Ls/s)
xs = [0.0; cumsum(fill(10s, L-1))]
mss = 0:0.05:2
Eps = map(mss) do ms
    model = Barriers.MIModel(2ms*s, xs, fill(locus, L))
    mean(model.Ep)
end
plot(mss, Eps)

