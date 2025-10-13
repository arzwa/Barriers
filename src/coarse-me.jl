"""
    CoarseModel

Coarse-grained window `mₑ` model.
"""
@with_kw struct CoarseModel{T,V}
    X :: Vector{Int}  # number of selected sites in window
    Δ :: Vector{T}    # winsizes in Morgan
    R :: Matrix{T} = winrecrates(Δ)  # between window recombination rates
    s :: V  # selection coefficient/vector of selection coefficients
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
    gs = map(1:length(X)) do i
        gi = map(1:length(X)) do j
            i == j ? log(gii(mod, i)) : log(gij(mod, i, j))
        end |> sum |> exp
    end  
end

_gets(s::Real, i) = s
_gets(s, i) = getindex(s, i)

function gii(model, i)
    @unpack Δ, m, u = model
    s = _gets(model.s, i)
    n = model.X[i]
    n == 0 && (return 1.0)
    δ  = Δ[i]/n
    xs = (δ/2):δ:(Δ[i]-δ/2)
    #diploids...
    #h = 0.5
    #f(x) = exp(-sum([s*h / (s*h + m + recrate(abs(x - xs[i]))) for i=1:n]))
    f(x) = exp(-sum([s / (s + m + recrate(abs(x - xs[i]))) for i=1:n]))
    g = quadgk(f, 0, Δ[i])[1]/Δ[i]
    #g * exp(-s*model.X[i])    # diploid model
    g
end

function _gii(model, i)
    @unpack Δ, m, u = model
    s = _gets(model.s, i)
    n = model.X[i]
    n == 0 && (return 1.0)
    δ  = Δ[i]/n
    xs = (δ/2):δ:(Δ[i]-δ/2)
    function f(x)
        r = findfirst(xi->xi > x, xs)
        r = isnothing(r) ? n+1 : r
        l = r - 1
        g = 1.0
        sleft = 0.0
        while l > 0
            rr = recrate(x - xs[l])
            g /= (1 + s/(rr + sleft))
            sleft += s
            l -= 1
        end
        srght = 0.0
        while r <= n
            rr = recrate(xs[r] - x)
            g /= (1 + s/(rr + srght))
            srght += s
            r += 1
        end
        return g
    end
    g = quadgk(f, 0, Δ[i])[1]/Δ[i]
    return g
end

gij(model, i, j) = exp(loggij(model, i, j))
function loggij(model, i, j)
    @unpack X, R, m, u = model
    s = _gets(model.s, j)
    #h = 0.5
    #-s*h*X[j]/(s*h + m + R[i,j]) - s*X[j]
    -s*X[j]/(s + m + R[i,j])   # haploid model
end

# for haploids...
function predict_divergence(model::CoarseModel, N, tol=1e-5)
    @unpack X, s, Δ, m, u = model
    s̄ = length(s) == 1 ? s[1] : mean(s)
    # get harmonic mean recombination rate
    r = harmonicmean_recrate(model)
    # predict divergence
    L = sum(X)
    gfun(Ep) = exp(-L*s̄*Ep/r)
    Ep = 1.0
    while true
        d = Wright(-2N*s̄, N*u, N*(m*gfun(Ep) + u), 0.5) 
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
