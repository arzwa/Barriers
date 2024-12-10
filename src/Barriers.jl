module Barriers

using Printf
using LinearAlgebra
using Sewall
using QuadGK
using AdaptiveProposals
using Distributions
using Serialization
using Parameters
using Accessors
using WrightDistribution
using Interpolations
using StatsBase
using Base.Iterators: partition
using SpecialFunctions

const Mb = 10^6
const kb = 10^3

include("geneticmap.jl")
include("windowed.jl")
include("aeschbacher_me.jl")
include("diffusion_me.jl")
include("strongbgs.jl")
#include("sample.jl")
include("sample3.jl")
include("sample4.jl")
include("sites_pr.jl")
include("infsite-likelihood.jl")

export Mb, kb
export WindowedChromosome, GeneticMap, AeschbacherModel, MIModel, MILocus
export me_profile, me


end # module Barriers
