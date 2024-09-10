using Plots, DataFrames, CSV, Distributions, Interpolations, Optim, Parameters
using AdaptiveProposals
using Sewall
plotsdefault()

df = CSV.read("/home/arthur_z/dev/Barriers/data/64bpX500blocks-rho-cds-sumstats-me-fit-range-comb.tsv", DataFrame)
df18 = filter(x->x[:chrom] == 18, df)

lmap = CSV.read("/home/arthur_z/dev/Barriers/data/linkage_map.cydno.cm.tsv", DataFrame)
lmap = filter(x->x[:Chromosome] == "chr18", lmap)

map18 = linear_interpolation(lmap[:,2], lmap[:,3])

plot(df18[!,:winMid], df18[!,:r_1000], legend=false, size=(900,200))
plot!(df18[!,:winMid], df18[!,:r_500], legend=false, size=(900,200))
plot!(df18[!,:winMid], df18[!,:r_250], legend=false, size=(900,200))

plot(df18[!,:winMid], df18[!,:f_st])
plot!(twinx(), df18[!,:winMid], df18[!,"m_cyd..mel"], color=2, legend=false, size=(400,200))

plot(df18[!,:winMid], df18[!,:r_1000], legend=false, size=(900,200))
plot!(twinx(), df18[!,:winMid], df18[!,:cds_250], color=2, legend=false)

data = df18[1:5:end,["winMid", "m_cyd..mel", "r_2000", "cds_250", "f_st", "N_mel"]]
rename!(data, "winMid"=>"x", "m_cyd..mel"=>"m", "r_2000"=>"r", "cds_250"=>"cds", "f_st"=>"fst") 
#data[:,:p] = [map18[x]/100 for x in data[:,:x]] 
data[:,"r"] ./= 100  # cM/Mb -> M/Mb

function getmap(data)
    between = [(data[i,:x] + data[i+1,:x])/2 for i=1:nrow(data)-1]
    xs = [0.0 ; between; data[end,:x]] ./ 1e6
    Wspan = [(xs[i], xs[i+1])  for i=1:length(xs)-1]
    ds = [(xs[i+1] - xs[i]) for i=1:length(xs)-1]  # length in Mb of each window
    gd = [ds[i] * data[i,:r] for i=1:length(ds)]  # Morgan
    gm = [0; cumsum(gd)]
    Mspan = [(gm[i], gm[i+1])  for i=1:length(gm)-1]
    gm = [(gm[i+1] + gm[i])/2 for i=1:length(gm)-1]
    ds, gd, gm, Mspan, Wspan
end

ds, gd, mpos, Mspan, Wspan = getmap(data)
data[:,:l] = ds    # window length in Mb
data[:,:d] = gd    # genetic length in M  
data[:,:p]  = mpos  # window midpoint in M
data[:,:Mspan] = Mspan  # window span in M
data[:,:Wspan] = Wspan  # window span in Mb

f = 1 ./ (1 .+ 11data[:,:N_mel] .* data[:,:m])
P1 = plot(data[:,:p], 1 ./ f, size=(500,300), 
    xlabel="map position", label="gIMble", ylabel="\$F_{ST}\$")
plot!(data[:,:p], 1 ./ data[:,:fst], label="observed")
P2 = scatter(1 ./ f, 1 ./ data[:,:fst], color=:black, ms=2, 
    xlabel="gIMble", ylabel="observed", legend=false) 
plot!(x->x, color=:gray, alpha=0.5)
plot(P1, P2, size=(600,250), margin=3Plots.mm)


plot(data[:,:x], data[:,:p]) 

plot(data[:,:p], data[:,:m], legend=false)
plot!(twinx(), data[:,:p], data[:,:cds], color=2, 
    size=(900,200), legend=false)

plot(data[:,:p], 1 ./ (1 .+ 4*1e6 .* data[:,:m]))
plot!(data[:,:p], data[:,:fst])

scatter(1 ./ (1 .+ 4*1e6 .* data[:,:m]), data[:,:fst], ms=2, color=:black,
    size=(400,350), legend=false, )

# build the genetic architecture for a given set of global parameters
function getarch(data, s, ν, μ, by=:l)  # by is he covariate by which we multiply nu to get the number of loci
    rs = map(eachrow(data)) do row
        nloci = ceil(Int64, ν*row[by])  # HACK: should do Poisson sampling?
        rij   = Sewall.haldane(row["d"]/nloci) 
        fill(rij, nloci)
    end |> x->vcat(x...)
    Architecture(fill(Locus(-s, 0.0, 0.0, μ, μ), length(rs)), rs[1:end-1])
end

# test
N = Ne2N(mean(df18[:,"N_mel"]),1)
m = mean(df18[:,"m_cyd..mel"])*2
ν = 15.0
s = 2e-4
μ = 1e-10
L = 500
A = getarch(data, s, ν, μ)
M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, HWLEMainland(ones(length(A))))
Ls = sum([-x.s1 for x in A])
p = mean.(equilibrium(M))
xs, ys = Sewall.meprofile(M, L)
# XXX should speed up calculation of the R matrix.
plot(xs, ys, legend=false, color=1, xlim=extrema(xs))
plot!(data[:,:p], data[:,:m], legend=false, size=(900,200))
vline!(Sewall.mappositions(A), color=:gray, lw=0.1)

P = plot(Sewall.mappositions(A), p, ylim=(0,1), size=(900,200), legend=false,
    color=1)

me_gimble = linear_interpolation(data[:,:p], data[:,:m])
m̂s = [me_gimble[x] for x in xs]
plot(xs, ys, legend=false, color=2)
plot!(twinx(), xs, m̂s, legend=false, size=(900,200))
plot!(twinx(), data[:,:p], data[:,:m], legend=false, color=3)

# likelihood function for the full data
_ℓhood(m, m̂, σ) = logpdf(Normal(m̂, σ), m)
ℓhood(ms, m̂s, σ) = sum(_ℓhood.(ms, m̂s, σ))

priors = (
    # a parameter of Gamma DFE
    # a=Exponential(2),
    # b parameter of Gamma DFE
    # b=Exponential(0.02/2),
    # selection coefficient
    s=Exponential(2e-4),
    # number of selected loci per CDS
    #ν=Exponential(1),
    ν=Exponential(10.0),
    # baseline migration rate
    m=Exponential(mean(df18[:,"m_cyd..mel"])*2),
    #m=Exponential(1e-6),
    μ=Dirac(1e-10),
    N=Dirac(Ne2N(mean(df18[:,"N_mel"]),1)),
    σ=Dirac(0.3)
)

logprior(priors, θ) = mapreduce(x->logpdf(x[1],x[2]), +, zip(priors, θ))

# To speed up: better initial allele frequencies for selected sites (mean of
# previous iteration?), more clever computation of R matrix
function mcmc_me(data, props, prior, n; every=10, L=1000)
    # initialize
    θ = map(rand, prior) 
    A = getarch(data, θ.s, θ.ν, θ.μ)
    M = MainlandIsland(Deme(N=θ.N, k=1, A=A), θ.m, 0.0, HWLEMainland(ones(length(A))))
    xs, ms = Sewall.meprofile(M, L)
    me_gimble = linear_interpolation(data[:,:d], data[:,:m])
    m̂s = [me_gimble[x] for x in xs] 
    pp = ℓhood(ms, m̂s, θ.σ) + logprior(priors, θ)
    pps = zeros(n)
    xxs = Vector{typeof(θ)}(undef, n)
    for k=1:n
        if k % every == 0 
#            @printf "%6d, %9.4f, %9.7f, %9.6f, %9.7f\n" k pp θ.s θ.ν θ.m
        end
        for sym in keys(props)
            @printf "%5s, %6d, %9.4f, %9.7f, %9.6f, %9.7f\n" sym k pp θ.s θ.ν θ.m
            prop = getfield(props, sym)
            x_, q = prop(getfield(θ, sym))
            θ_ = reconstruct(θ, sym=>x_)
            A_ = getarch(data, θ_.s, θ_.ν, θ_.μ)
            y_ = HWLEMainland(ones(length(A_)))
            M_ = MainlandIsland(Deme(N=θ_.N, k=1, A=A_), θ_.m, 0.0, y_)
            pp_, ms_ = try
                _, ms_ = Sewall.meprofile(M_, L)
                # XXX `getarch` computes `R`, but this is not needed, since we
                # compute it again but with neutral loci in `meprofile`.
                # XXX `meprofile` is not yet optimized.
                pp_ = ℓhood(ms_, m̂s, θ_.σ) + logprior(priors, θ_)
                pp_, ms_
            catch
                @warn "Numerical error"
                pp_ = -Inf
                pp_, similar(ms)
            end
            if log(rand()) < pp_ - pp + q
                pp = pp_
                θ = θ_
                accept!(prop)
            end
        end
        pps[k] = pp
        xxs[k] = θ
    end
    return pps, xxs
end

proposals = (;[sym=>PositiveProposal(kernel=Normal(0,1)) for sym in [:s, :m]]...)

pps, xxs = mcmc_me(data, proposals, priors, 1000, L=100, every=1)

# 1. We have global parameters, with a global contribution to the prior (dfe
# parameters, number of selected loci, baseline m). In the simplest model we
# assume a fixed selection coefficient. So we have `s`, `m` and `ν`
# 2. We have window-specific parameters (selection coefficients)
# 3. Resampling a global parameter changes all the window-specific
# contributions
# 4. But this is also the case for the window specific parameters -> changes
# the equilibrium frequencies everywhere, and hence the gff evereywhere, not
# only in that particular window (note though that in an MCMC algorithm one
# could initialize the fixed point iteration with the previous allele
# frequencies, which should reduce at least the computational cost of
# calculating the gffs. 
# 5. One could start by assuming equally spaced selected loci

function optimize_me(data, init; μ=1e-10, N=1e6, L=100)
    # a bit silly, but its just to get the `xs`
    s, ν, m = init
    A = getarch(data, s, ν, μ)
    y = HWLEMainland(ones(length(A)))
    M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
    xs, ms = Sewall.meprofile(M, L)
    me_gimble = linear_interpolation(data[:,:p], data[:,:m])
    m̂s = [me_gimble[x] for x in xs] 
    function obj(x) 
        s, ν, m = exp.(x)
        A = getarch(data, s, ν, μ)
        y = HWLEMainland(ones(length(A)))
        M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
        try 
            _, ms = Sewall.meprofile(M, L)
            ssd = sum(((ms .- m̂s) ./ m) .^ 2)
            @info x, ssd
            ssd
        catch
            Inf
        end
    end
    optimize(obj, log.(init), NelderMead())
end

x = map(rand, priors)
init = [x.s, x.ν, x.m]
res = optimize_me(data, init, L=200)

# res.minimizer = [-16.480947909617328, 2.1333585066871206, -13.851107328686453]
# exp.(res.minimizer) = [6.956892079633877e-8, 8.44317570909929, 9.65029343124853e-7]
x̂ = exp.(res.minimizer)
A = getarch(data, x̂[1], x̂[2], 1e-10)
y = HWLEMainland(ones(length(A)))
M = MainlandIsland(Deme(N=9e5, k=1, A=A), x̂[3], 0.0, y)
xs, ms = Sewall.meprofile(M, L)
me_gimble = linear_interpolation(data[:,:p], data[:,:m])

plot(xs, ms, size=(900,200))
plot!(xs, [me_gimble[x] for x in xs])
# so this just infers no selection but a migration rate which is E[m_e]

# fix s
function optimize_me_2(data, init; s=2e-4, μ=1e-10, N=1e6, L=500, 
        kwargs...)
    # a bit silly, but its just to get the `xs`
    ν, m = init
    A = getarch(data, s, ν, μ)
    y = HWLEMainland(ones(length(A)))
    M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
    xs, ms = Sewall.meprofile(M, L)
    me_gimble = linear_interpolation(data[:,:p], data[:,:m])
    m̂s = [me_gimble[x] for x in xs] 
    function obj(x) 
        ν, m = exp.(x)
        A = getarch(data, s, ν, μ)
        y = HWLEMainland(ones(length(A)))
        M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
        try 
            _, ms = Sewall.meprofile(M, L)
            ssd = sum(((ms .- m̂s) ./ m) .^ 2)
            @info x, ssd
            ssd
        catch
            @warn "Numerics"
            Inf
        end
    end
    optimize(obj, log.(init), NelderMead(),
        Optim.Options(; kwargs...))
end

ress = map([20/1e6, 10/1e6, 2/1e6]) do s
    init = exp.([2.0, -12.0])
    optimize_me_2(data, init, N=1e6, μ=1e-10, s=s, L=200, iterations=200)
end
    
init = exp.([2.0, -12.0])
optimize_me_2(data, init, N=1e6, μ=1e-10, s=10.0/1e6, 
    L=100, iterations=200)

#x̂ = exp.(res2.minimizer)
x̂ = exp.([4.053, -13.812])
A = getarch(data, s, x̂[1], 1e-10)
y = HWLEMainland(ones(length(A)))
M = MainlandIsland(Deme(N=9e5, k=1, A=A), x̂[2], 0.0, y)
xs, ms = Sewall.meprofile(M, L)
me_gimble = linear_interpolation(data[:,:p], data[:,:m])

plot(xs, ms, size=(900,200))
plot!(xs, [me_gimble[x] for x in xs])


# fix ν
function optimize_me_3(data, init; ν=2.0, μ=1e-10, N=1e6, L=500, kwargs...)
    # a bit silly, but its just to get the `xs`
    s, m = init
    A = getarch(data, s, ν, μ, :cds)
    y = HWLEMainland(ones(length(A)))
    M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
    xs, ms = Sewall.meprofile(M, L)
    me_gimble = linear_interpolation(data[:,:p], data[:,:m])
    m̂s = [me_gimble[x] for x in xs] 
    function obj(x) 
        ν, m = exp.(x)
        A = getarch(data, s, ν, μ)
        y = HWLEMainland(ones(length(A)))
        M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
        try 
            _, ms = Sewall.meprofile(M, L)
            ssd = sum(((ms .- m̂s) ./ m) .^ 2)
            @info x, ssd
            ssd
        catch
            @warn "Numerics"
            Inf
        end
    end
    optimize(obj, log.(init), NelderMead(),
        Optim.Options(; kwargs...))
end

ν = 1.0
N = 1e6
init = [1/N, 0.2/N]
res3 = optimize_me_3(data, init, ν=ν, N=N, μ=1e-10, L=100, iterations=200)
x̂ = exp.(res3.minimizer)
@info x̂[1]*N

A = getarch(data, x̂[1], ν, 1e-10)
y = HWLEMainland(ones(length(A)))
M = MainlandIsland(Deme(N=N, k=1, A=A), x̂[2], 0.0, y)
xs, ms = Sewall.meprofile(M, L)
me_gimble = linear_interpolation(data[:,:p], data[:,:m])
plot(xs, ms, size=(900,200))
plot!(xs, [me_gimble[x] for x in xs])


# Sampler, again. For fixed ν, infer s, m and sample locations.

function sample_window(data, ν, i; α=10.0)
    # sample number of loci and their positions for window i
    winlen  = data[i,:l]      # length in Mb
    x1, _   = data[i,:Mspan]  # coordinates of window
    windist = data[i,:d]      # length of window in M
    dd = Poisson(ν*winlen)
    nloci = rand(dd)
    ℓ = logpdf(dd, nloci)
    if nloci > 0 
        dp = Dirichlet(nloci, α)
        pos01 = rand(dp)
        pos = x1 .+ cumsum(pos01) .* windist 
        ℓ += logpdf(dp, pos01)
    else
        pos = Float64[]
    end
    return pos, ℓ
end

function sample_arch(data, ν; α=10.0)
    map(i->sample_window(data, ν, i, α=α), 1:nrow(data))
end

function resample_window(windows, i, data, ν; α=10.0)
    prev = windows[i]
    pos, l = sample_window(data, ν, i, α=α)
    windows[i] = (pos, l)
    return windows, prev
end
    
function recrates(windows)
    pos = vcat(first.(windows)...)
    Sewall.ratematrix(pos)
end


function mprediction(windows, θ, nn)
    @unpack ν, s, m, σ, μ, N = θ
    R = recrates(windows)
    L = size(R,1)
    A = Architecture([Locus(-s, 0.0, 0.0, μ, μ) for i=1:L], R)
    y = HWLEMainland(ones(L))
    M = MainlandIsland(Deme(N=N, k=1, A=A), m, 0.0, y)
    return Sewall.meprofile(M, nn)
end

getfst(me, N) = 1 ./ (1 .+ 4N .* me) 

logprior(priors, θ) = mapreduce(x->logpdf(x[1],x[2]), +, zip(priors, θ))

# likelihood function for the full data
_ℓhood(m, m̂, σ) = logpdf(Normal(m̂, σ), m)
ℓhood(ms, m̂s, σ) = sum(_ℓhood.(ms, m̂s, σ))

function sampler(data, θ, prior, props, n; α=10.0, nn=100, nw=10)
    @unpack ν, s, m, σ, μ, N = θ
    me_gimble = linear_interpolation(data[:,:p], data[:,:m])
    windows = sample_arch(data, ν; α=α)
    xs, ms = mprediction(windows, θ, nn)
    m̂s = [me_gimble[x] for x in xs] 
    pp = ℓhood(ms, m̂s, θ.σ) + logprior(prior, θ) + sum(last.(windows))
    pps = zeros(n)
    xxs = Vector{typeof(windows)}(undef, n)
    θs  = Vector{typeof(θ)}(undef, n)
    xxs[1] = copy(windows)     
    pps[1] = pp
    for k=2:n
        @printf "(%4d) π=%12.8f Ns=%9.5f Nm=%9.5f\n" k pp θ.N*θ.s θ.N*θ.m
        for sym in keys(props)
            prop = getfield(props, sym)
            x_, q = prop(getfield(θ, sym))
            θ_ = reconstruct(θ, sym=>x_)
            windows_ = sample_arch(data, θ_.ν; α=α)
            _, ms = mprediction(windows_, θ_, nn)
            pp_ = ℓhood(ms, m̂s, θ_.σ) + logprior(prior, θ_) + sum(last.(windows_)) + q
            if log(rand()) < pp_ - pp
                θ = θ_
                windows = windows_
                pp = pp_
                accept!(prop)
            end
        end
        @unpack ν, s, m, σ, μ, N = θ
        for i=rand(1:nrow(data), nw)
            #@printf "-> %12.8f\n" pp
            windows, prev = resample_window(windows, i, data, ν; α=α)
            _, ms = mprediction(windows, θ, nn)
            pp_ = ℓhood(ms, m̂s, θ.σ) + logprior(prior, θ) + sum(last.(windows))
            if log(rand()) < pp_ - pp
                pp = pp_
            else
                windows[i] = prev
            end
        end
        θs[k]  = θ 
        pps[k] = pp
        xxs[k] = copy(windows)
    end
    return θs, pps, xxs
end

priors = (
    ν=Dirac(10.0),
    s=Exponential(1e-5), 
    m=Exponential(1e-6),
    σ=Dirac(1e-5),
    μ=Dirac(1e-10),
    N=Dirac(1e6)
)
x0 = map(rand, priors)
props = (
    s=PositiveProposal(kernel=Normal(0,0.1)), 
    m=PositiveProposal(kernel=Normal(0,0.1)))
sampler(data, x0, priors, props, 100)

function sampler_fst(data, θ, prior, props, n; α=10.0, nn=100, nw=10)
    @unpack ν, s, m, σ, μ, N = θ
    fst = linear_interpolation(data[:,:p], data[:,:fst])
    windows = sample_arch(data, ν; α=α)
    xs, ms = mprediction(windows, θ, nn)
    Yobs  = [fst[x] for x in xs] 
    Ypred = getfst(ms, N)
    obj(y, y_) = sum(y .- y_)
    @info obj(Yobs, Ypred)
    pp = logpdf(Normal(0, σ), obj(Yobs, Ypred)) + 
        logprior(priors, θ) + sum(last.(windows))
    pps = zeros(n)
    xxs = Vector{typeof(windows)}(undef, n)
    θs  = Vector{typeof(θ)}(undef, n)
    xxs[1] = copy(windows)     
    pps[1] = pp
    for k=2:n
        @printf "(%4d) π=%12.8f Ns=%9.5f Nm=%9.5f\n" k pp θ.N*θ.s θ.N*θ.m
        for sym in keys(props)
            prop = getfield(props, sym)
            x_, q = prop(getfield(θ, sym))
            θ_ = reconstruct(θ, sym=>x_)
            windows_ = sample_arch(data, θ_.ν; α=α)
            _, ms = mprediction(windows_, θ_, nn)
            Ypred = getfst(ms, N)
            pp_ = logpdf(Normal(0, σ), obj(Yobs, Ypred)) + 
                logprior(priors, θ_) + sum(last.(windows_)) + q
            if log(rand()) < pp_ - pp
                θ = θ_
                windows = windows_
                pp = pp_
                accept!(prop)
            end
        end
        @unpack ν, s, m, σ, μ, N = θ
        for i=rand(1:nrow(data), nw)
            #@printf "-> %12.8f\n" pp
            windows, prev = resample_window(windows, i, data, ν; α=α)
            _, ms = mprediction(windows, θ, nn)
            Ypred = getfst(ms, N)
            pp_ = logpdf(Normal(0, σ), obj(Yobs, Ypred)) +
                logprior(prior, θ) + sum(last.(windows))
            if log(rand()) < pp_ - pp
                pp = pp_
            else
                windows[i] = prev
            end
        end
        θs[k]  = θ 
        pps[k] = pp
        xxs[k] = copy(windows)
    end
    return θs, pps, xxs
end

priors = (
    ν=Dirac(5.0),
#    s=Dirac(1e-5), 
    s=Exponential(1e-5),
#    m=Exponential(1e-6),
    m=Dirac(0.5/1e6),
    σ=Dirac(1.0),
    μ=Dirac(1e-10),
    N=Dirac(1e6)
)
x0 = map(rand, priors)
props = (
    s=PositiveProposal(kernel=Normal(0,1)),)
#    m=PositiveProposal(kernel=Normal(0,1)),)
ts, ps, xs = sampler_fst(data, x0, priors, props, 1000, nw=20)

plot(data[:,:p], data[:,:fst])
for i=2:10:100
    x, y = mprediction(xs[i], ts[i], 100)
    plot!(x, getfst(y, ts[i].N))
end
plot!()

hcat(map(i->map(x->length(x[i][1]), xs), 1:nrow(data))...)

plot(data[:,:p], data[:,:fst])
plot!(twinx(), data[:,:p], vec(mean(hcat(map(i->map(x->length(x[i][1]), xs),
    1:nrow(data))...), dims=1)), color=2, legend=false)


# Either we get a proper data model, or we implement a proper SA-like algorithm
# for fitting...
