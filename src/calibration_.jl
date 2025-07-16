module Calibration

using Barriers
using PyCall, StatsBase, Distributions, Parameters
using ProgressMeter, Combinatorics, Optim, QuadGK

# parameters are {NA, NB, u, m, r}
struct MigrationCalibration{T,V}
    prior :: V 
    NA :: T
    NB :: T
    u  :: T
    r  :: T    # per bp
    L  :: Int  # bp
    na :: Int
    nb :: Int
end

function simulations(C::MigrationCalibration, nrep)
    @unpack NA, NB, u, r, L, na, nb = C
    @showprogress desc="Simulating $nrep replicates" map(1:nrep) do _
        m = rand(C.prior)
        M = TwoDemeCoalescent(m=m, NA=NA, NB=NB, u=u, r=r, L=L, na=na, nb=nb)
        y = randvars(M)
    end
end

function simulation()

end
