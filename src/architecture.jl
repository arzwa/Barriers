abstract type Locus end

struct Architecture{L<:Locus,T<:Real}
    loci :: Vector{L}
    xs   :: Vector{T}  # map positions
    R    :: Matrix{T}  # pairwise recombination rates
end

Architecture(loci, xs::AbstractVector) = Architecture(loci, xs, recrates(xs))

Base.length(A::Architecture) = length(A.loci)
Base.getindex(A::Architecture, i) = A.loci[i]

struct DiploidLocus{T<:Real} <: Locus
    s :: T  # selection coefficient 
    h :: T  # dominance coefficient
    u :: T  # mutation rate
end

HaploidLocus(s, u) = DiploidLocus(2s, 0.5, u)
