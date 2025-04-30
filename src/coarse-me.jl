"""
    CoarseModel

Coarse-grained window `mₑ` model.
"""
@with_kw struct CoarseModel{T}
    X :: Vector{Int}  # number of selected sites in window
    Δ :: Vector{T}    # winsizes in Morgan
    R :: Matrix{T} = winrecrates(Δ)  # between window recombination rates
    s :: T  # selection coefficient
    m :: T  # migration rate
    u :: T  # mutation rate
    λ :: T  # coal. rate (inverse pop size)
end

# XXX mutation

# accounting for partial divergence in a crude way (using harmonic mean
# recombination rate to calculate a single Ep across the barrier). 
function gff(model::CoarseModel)
    @unpack X, Δ, R, s, m, u, λ = model 
    mod = if λ == 0.  
        # complete divergence == no drift == coalescence rate 0
        model
    else
        # account for partial divergence (drift)
        div = predict_divergence(model, 1/λ)
        CoarseModel(X, Δ, R, s*div, m, u, λ)
    end
    map(1:length(X)) do i
        gi = map(1:length(X)) do j
            i == j ? log(gii(mod, i)) : log(gij(mod, i, j))
        end |> sum |> exp
    end  
end

function gii(model, i)
    @unpack s, Δ, m, u = model
    n = model.X[i]
    n == 0 && (return 1.0)
    x = range(0, Δ[i], n+1)
    xs = [(x[i] + x[i+1])/2 for i=1:n]
    f(x) = exp(-sum([s / (s + m + recrate(abs(x - xs[i]))) for i=1:n]))
    g = quadgk(f, 0, Δ[i])[1]/Δ[i]
end

gij(model, i, j) = exp(loggij(model, i, j))
function loggij(model, i, j)
    @unpack X, s, R, m, u = model
    -s*X[j]/(s + m + R[i,j])
end

# for haploids...
function predict_divergence(model::CoarseModel, N, tol=1e-5)
    @unpack X, s, Δ, m, u = model
    # get harmonic mean recombination rate
    r = harmonicmean_recrate(model)
    # predict divergence
    L = sum(X)
    gfun(Ep) = exp(-L*s*Ep/r)
    Ep = 1.0
    while true
        d = Wright(-2N*s, N*u, N*(m*gfun(Ep) + u), 0.5) 
        Ep_ = mean(d)
        abs(Ep_ - Ep) < tol && return Ep_
        Ep = Ep_
    end
end

# These are the implied map positions, assuming equally spaced loci within
# each window
function mappositions(model::CoarseModel)
    @unpack X, Δ = model
    xs = Float64[]
    x = 0.0 
    for i=1:length(X)
        # X[i] is the number of selected sites in the window
        # Δ[i] is the window size in Morgans
        dx = Δ[i]/(X[i]+1)
        x += dx
        for j=1:X[i]
            push!(xs, x)
            x += dx
        end
    end
    return xs
end

function harmonicmean_recrate(model::CoarseModel)
    xs = mappositions(model)
    r  = 0.0
    for i=1:length(xs)
        ri = 0.0
        for j=1:length(xs)
            j == i && continue
            ri += 1 / recrate(abs(xs[i] - xs[j]))
        end
        r += 1 / (ri / (length(xs)-1))
    end
    r / length(xs)
end

# Calculate recombination rates between window midpoints
function winrecrates(Δ::Vector)
    n = length(Δ)
    R = zeros(n,n)
    for i=1:n-1
        for j=i+1:n
            d = Δ[i]/2 + sum(Δ[i+1:j-1]) + Δ[j]/2
            R[i,j] = R[j,i] = recrate(d)
        end
    end
    return R
end
