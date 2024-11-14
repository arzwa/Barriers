using CSV, DataFrames
using Barriers

df = CSV.read("/home/arthur_z/dev/Barriers/data/64bpX500blocks-rho-cds-sumstats-me-fit-range-comb.tsv", DataFrame)
df18 = filter(x->x[:chrom] == 18, df)

lmap = CSV.read("/home/arthur_z/dev/Barriers/data/linkage_map.cydno.cm.tsv", DataFrame)
lmap = filter(x->x[:Chromosome] == "chr18", lmap)
lmap = lmap[.!nonunique(lmap[:,[:cM]]),:]
gmap = GeneticMap(collect(zip(lmap[:,:End] .+ 1, lmap[:,:cM] ./ 100)))

# take this as the end of the chromosome (a window is approximately 100kb) 
chr_end = maximum(df18[:,:winMid])+50kb

# issue!
chr_end > Barriers.physlength(gmap)

# we should extrapolate based on the average recombination rate?
gmap = Barriers.extrapolate_to(gmap, chr_end)
chr_end == Barriers.physlength(gmap)

# nonoverlapping windows
data = df18[1:5:end,:]
xs = data[:,:winMid]
zs = [1 ; floor.(Int64,[(xs[i] + xs[i-1])/2 for i=2:length(xs)]) ; chr_end+1]
wins = [zs[i]:zs[i+1]-1 for i=1:length(zs)-1]

# the object
d = WindowedChromosome(gmap, wins)
meobs = [win => m for (win, m) in zip(first.(d.data), data[:,"m_cyd..mel"])]
