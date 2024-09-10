module Barriers

using Printf
using LinearAlgebra
using Sewall
using QuadGK
using AdaptiveProposals
using Distributions
using Serialization
using Parameters
using WrightDistribution
using Interpolations
using Base.Iterators: partition

const Mb = 10^6
const kb = 10^3

include("windows.jl")
include("me.jl")

export Mb, kb


end # module Barriers
