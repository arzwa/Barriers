abstract type Locus end

struct Architecture{L<:Locus,T<:Real}
    loci :: Vector{L}
    xs   :: Vector{T}  # map positions
    R    :: Matrix{T}  # pairwise recombination rates
end

Base.length(A::Architecture) = length(A.loci)
Base.getindex(A::Architecture, i) = A.loci[i]

struct DiploidLocus{T<:Real} <: Locus
    s :: T  # selection coefficient 
    h :: T  # dominance coefficient
    u :: T  # mutation rate
end

