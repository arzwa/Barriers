using Barriers, Distributions

ν = 1e-5
m = 1e-4
s = 1e-3
σ = 1
n = length(d)
X = [rand(Poisson(ν*length(d[i][1]))) for i=1:n]
sum(X)*s

state  = Barriers.State((m=m, s=s, ν=ν, σ=σ), X, Vector{Float64}(undef, n))
mepred = [Barriers.predict_me(state, d, i) for i=1:n]
state  = reconstruct(state, mepred=mepred)

xx = [last(d[i][1]) for i=1:n]
plot(xx, mepred ./ m, size=(800,200), line=:steppost)
plot!(twinx(), xx, X, color=:gray, fill=true, alpha=0.2, line=:steppost)

# Given all other parameters, resample number of selected sites in the window
# without doing an MH step. f(X|m) ∝ f(m|X)f(X)
XX = [Barriers.probs(state, mepred[i], d, i) for i=1:n]
X_ = map(XX) do p
    sample(0:length(p)-1, Weights(p))
end
sum(X_)  # makes sense

plot(xx, mepred ./ m, size=(800,200), line=:steppost)
plot!(twinx(), xx, X, color=:gray, fill=true, alpha=0.2, line=:steppost)
plot!(twinx(), xx, X_, color=:red, fill=true, alpha=0.2, line=:steppost)

# Now should implement full inference under this approximate model, i.e.
# estimate ν, s, m for given σ and a bunch of priors.

# Then should implement the more advanced version, for the same coarse
# architecture: we can calculate an m_e, so we can put it in a fixed-point
# iteration?
