"""
    AeschbacherModel(m, s::Vector, xs::Vector)

- `m`  migration rate
- `s`  vector of selection coefficients (haploid selection)
- `xs` vector of map positions (in Morgan)
"""
struct AeschbacherModel{T,V}
    m :: T          # migration rate
    s :: V          # selection coefficients
    xs:: Vector{T}  # map positions of the selected loci
    n :: Int64      # consider n MSPs on either side
    function AeschbacherModel(m::T, s::V, xs::Vector{T}; n=100) where {T,V} 
        new{T,V}(m, s, sort(xs), n)
    end
end

gets(M::AeschbacherModel{T,V}, i) where {T,V<:AbstractVector} = M.s[i]
gets(M::AeschbacherModel{T,V}, _) where {T,V<:Real} = M.s

# calculate me at **map** position x
me(model::AeschbacherModel, x) = model.m*gff(model, x)
function gff(model::AeschbacherModel, x)
    @unpack m, xs = model
    loggff = 0.0
    function _recursel(a, i)
        (i == length(xs) + 1 || xs[i] > x) && return 0.0
        b = _recursel(a, i+1)
        r = recrate(abs(xs[i] - x))
        loggff -= log(1 - gets(model, i)/(r + b))
        return a - b  # Aeschbacher expressions are for positive s...
    end
    function _recurser(a, i)
        (i == 0 || xs[i] < x) && return 0.0
        b = _recurser(a, i-1)
        r = recrate(abs(xs[i] - x))
        loggff -= log(1 - gets(model, i)/(r + b))
        return a - b
    end
    left = findlast(z->z < x, xs)
    left = isnothing(left) ? 1 : left
    xl = max(1, left - model.n + 1)
    xr = min(length(xs), left + model.n)
    _recursel(0.0, xl)
    _recurser(0.0, xr)
    return exp(loggff)
end

