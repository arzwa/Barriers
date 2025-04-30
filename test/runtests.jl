using Test, Barriers

@testset "Haploid multilocus diffusion" begin
    # Compare against Sachdeva fig 1A.
    L = 80 
    s = 0.02
    u = s*0.005
    N = 100
    # haploid selection sₑ = 2s!
    loci = fill(Barriers.DiploidLocus(2s, 0.5, u), L)
    A  = Barriers.Architecture(loci, fill(Inf, L), fill(0.5, L, L))
    M = Barriers.MainlandIslandModel(arch=A, m=0.5*s, N=N)
    p = Barriers.Equilibrium(M).Ep[1]
    ps = map(0.0:0.1:0.8) do ms
        M = Barriers.MainlandIslandModel(arch=A, m=ms*s, N=N)
        p = Barriers.Equilibrium(M).Ep[1]
        (ms, p)
    end
    # something weird happening with mutation rates
    @test 0.75 < ps[4] < 0.85
    @test 0.03 < ps[5] < 0.13
end

@testset "Diploid diffusion" begin
    # NOTE: when comparing against Zwaenepoel 2024, note that there Ne*s uses
    # the haplodiplontic Ne, which is the effective number of gene
    # copies/haplotypes. i.e. Ne*s is 2N*s, so when they used Ne*s=20, this
    # means N = 20/2s here (and recall that in the diffusion one should use
    # then 2N to account for diploidy... factors of two everywhere...)
    Ls = 2.0
    h  = 0.0
    Ns = 20.
    s  = 0.01
    L  = ceil(Int64, Ls/s)
    N  = ceil(Int64, Ns/2s)
    u  = s/100
    ms = 0.05:0.02:1.25
    loci = fill(Barriers.DiploidLocus(s, h, u), L)
    A  = Barriers.Architecture(loci, fill(Inf, L), fill(0.5, L, L))
    ps = map(ms) do m
        M = Barriers.MainlandIslandModel(arch=A, m=m*s, N=2N)
        p = Barriers.Equilibrium(M).Ep[1]
    end
end

@testset "Diploid diffusion, heterogeneous" begin
    sh = [(0.01784812271919177, 0.4774294948245706, 4.0e-5), (0.0033427119330979643, 0.570189199625785, 4.0e-5), (0.010607115649965396, 0.6310222868077683, 4.0e-5), (0.022763822160632906, 0.48458180244607213, 4.0e-5), (0.05049476550302885, 0.4480014457384327, 4.0e-5), (0.03194435866629766, 0.4756731822042609, 4.0e-5), (0.07330305884944699, 0.6529935597867397, 4.0e-5), (3.210000999384292e-5, 0.5310636105775212, 4.0e-5), (0.005868388243948209, 0.6095911499191312, 4.0e-5), (0.016424374479843542, 0.525968689252523, 4.0e-5), (0.08964093084228708, 0.516122424198271, 4.0e-5), (0.02075649636886527, 0.46741728899144447, 4.0e-5), (0.05067260306401162, 0.4989244039253686, 4.0e-5), (0.04069671843190063, 0.5398572298180582, 4.0e-5), (0.0012214534848979286, 0.35841515128742274, 4.0e-5), (0.027925029749604955, 0.4979517510664406, 4.0e-5), (0.030130384377393397, 0.612632743018635, 4.0e-5), (0.02093217543077693, 0.6437673667024638, 4.0e-5), (0.00914103450376985, 0.4382767555247339, 4.0e-5), (0.00655628446507645, 0.5008992856640958, 4.0e-5), (0.046736867174036655, 0.5186198422685501, 4.0e-5), (0.05331576807561197, 0.6598242015139981, 4.0e-5), (0.12398600790745126, 0.4643787501034338, 4.0e-5), (0.01946591833674737, 0.5737286074436032, 4.0e-5), (0.026076550512897145, 0.5252767126184434, 4.0e-5), (0.025232799892400806, 0.4984936556350982, 4.0e-5), (0.01415806576677534, 0.5759593798144186, 4.0e-5), (0.08325950323809264, 0.498276635441459, 4.0e-5), (0.12067295704266764, 0.5372682578960544, 4.0e-5), (0.01698202751117508, 0.3354807377823536, 4.0e-5), (0.03275003065794742, 0.513186460402217, 4.0e-5), (0.007600866379694988, 0.48186638522838204, 4.0e-5), (0.034299837243203256, 0.4691255042815275, 4.0e-5), (0.009888342559392376, 0.35236567456042456, 4.0e-5), (0.015216181735916086, 0.5017231656973983, 4.0e-5), (0.001261667488499526, 0.4833638536004813, 4.0e-5), (0.0035093178116451813, 0.5206392051005292, 4.0e-5), (0.004099453813979295, 0.4981380422745971, 4.0e-5), (0.048019685662431935, 0.5183163399161754, 4.0e-5), (0.025732978608533526, 0.6728989979276748, 4.0e-5), (0.0054882490118892055, 0.4935469976184749, 4.0e-5), (0.06166203348959155, 0.3683059275897485, 4.0e-5), (0.0009524138006460401, 0.6377752054768278, 4.0e-5), (0.03246672589135768, 0.7805083033485235, 4.0e-5), (0.011033689894898802, 0.7841746335292585, 4.0e-5), (0.05189759779340604, 0.5138957656384034, 4.0e-5), (0.010290353639396628, 0.5438048461968003, 4.0e-5), (0.0010861435670569393, 0.5089588480281182, 4.0e-5), (0.10921357517844482, 0.48868550039515385, 4.0e-5), (0.027833951835916006, 0.4613615311132817, 4.0e-5)]
    loci = [Barriers.DiploidLocus(x...) for x in sh] 
    N = 500
    L = length(loci)
    A = Barriers.Architecture(loci, fill(Inf, L), fill(0.5, L, L))
    M = Barriers.MainlandIslandModel(arch=A, m=0.4*0.02, N=N)
    eq = Barriers.Equilibrium(M)
    @test all(eq.Ep[1:5] .≈ [0.5101708662142199, 0.027031413297720276, 0.10976065183051234, 0.7094968247768757, 0.8879364925871525])
end

# debug...
L = 80 
s = 0.02
u = s*0.005
N = 100
loci = fill(Barriers.DiploidLocus(2s, 0.5, u), L)
A = Barriers.Architecture(loci, fill(Inf, L), fill(0.5, L, L))
M = Barriers.MainlandIslandModel(arch=A, m=0.50*s, N=N)
p = Barriers.Equilibrium(M).Ep[1]

#L  = 40
#Ns = 4.0
#Ls = 0.8
#s  = 0.02
#N  = Ns/s
#u  = s*0.005
#loci = fill(Barriers.DiploidLocus(2s, 0.5, u), L)
#A  = Barriers.Architecture(loci, fill(Inf, L), fill(0.5, L, L))
#res = map(0:0.05:1.0) do ms
#    M = Barriers.MainlandIslandModel(arch=A, m=ms*s, N=N)
#    p = Barriers.Equilibrium(M).Ep[1]
#    ms, p
#end
#plot(res)







