@with_kw struct BPTwoPop{T,V,R<:RecombinationMap,N}
    m12 :: T
    m21 :: T
    s1  :: Vector{T}
    s2  :: Vector{T}
    xs  :: Vector{V}     # map positions of the selected loci
    Ne1 :: N 
    Ne2 :: N 
    u   :: T
    rm  :: R
end

function Equilibrium(model::BPTwoPop; 
        Ep1=ones(length(model.xs)), 
        Ep2=ones(length(model.xs)), 
        tol=1e-6, nmax=Inf, α=1.)
    @unpack rm, xs, s1, s2, m12, m21, Ne1, Ne2, u = model
    it = 0
    #ps = [[Ep1 ; Ep2]]
    while true  # for small m ps2 ≈ 1, ps1 ≈ 0
        div = Ep2 .- (1 .- Ep1)
        gs1 = gff(rm, xs, s1 .* div)
        gs2 = gff(rm, xs, s2 .* div)  
        _Ep = map(1:length(xs)) do i
            d1 = Wright(-2Ne1*s1[i], 
                Ne1*(m21*gs1[i]*(1-Ep2[i]) + u), 
                Ne1*(m21*gs1[i]*Ep2[i] + u))
            p1 = mean(d1)
            d2 = Wright(-2Ne2*s2[i], 
                Ne2*(m12*gs2[i]*(1-Ep1[i]) + u), 
                Ne2*(m12*gs2[i]*Ep1[i] + u))
            p2 = mean(d2)
            p1, p2
        end
        _Ep1 = first.(_Ep)
        _Ep2 = last.(_Ep)
        _Ep1 = (Ep1 .* (1-α)) .+ (α .* _Ep1)
        _Ep2 = (Ep2 .* (1-α)) .+ (α .* _Ep2)
        if norm(Ep1 .- _Ep1) + norm(Ep2 .- _Ep2) < tol || it > nmax
            #return permutedims(hcat(ps...))
            return Equilibrium(model, [Ep1 Ep2], fill(NaN, 2, length(Ep1)))
        end
        Ep1 = _Ep1
        Ep2 = _Ep2
        #push!(ps, [Ep1 ; Ep2])
        it += 1
    end
end

function gff(model::Equilibrium{M}, x) where M<:BPTwoPop
    @unpack m12, m21, s1, s2, xs, rm = model.model
    div = model.Ep[:,2] .- (1 .- model.Ep[:,1])
    g1 = gff(rm, x, xs, div .* s1)
    g2 = gff(rm, x, xs, div .* s2)
    g1, g2
end

function me(model::Equilibrium{M}, x) where M<:BPTwoPop
    @unpack m12, m21 = model.model
    g1, g2 = gff(model, x)
    g1 * m21, g2 * m12
end

