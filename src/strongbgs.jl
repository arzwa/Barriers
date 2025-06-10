# Hudson & Kaplan like 'B' maps 
# Haploid deleterious locus
struct DelLocus{T}
    s :: T
    u :: T
end

struct BGSModel{T,V}
    xs   :: Vector{T}  # positions of sites under purifying selection
    loci :: Vector{DelLocus{V}}
end

function B(model::BGSModel, x)
    @unpack xs, loci = model
    lB = 0.0
    for j=1:length(loci)
        lB += locuseffect(loci[j], recrate(abs(x - xs[j])))
    end
    return exp(-lB)
end

function locuseffect(l::DelLocus, rij)
    @unpack u, s = l
    u / (s*(1 + rij*(1-s)/s)^2)
end

