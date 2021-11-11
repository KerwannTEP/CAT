include("Args.jl") # Parsing the command-line
##################################################
# General parameters
##################################################
"""
    PARALLEL

Determining if the code is run in parallel.
"""
const PARALLEL = parsed_args["parallel"]
if ((PARALLEL != "yes") && (PARALLEL != "no"))
    error("ERROR: UNKNOWN PARALLEL") # Unknown parallel procedure
end


"""
    nbK_default

Default number of sampling points for the 3D integrals.
"""
const nbK_default = parsed_args["nbK"]

"""
    nbAvr_default

Default number of sampling points for the orbit-averaging integrals.
"""
const nbAvr_default = parsed_args["nbAvr"]

"""
    qCalc

Anisotropy parameter q for the Plummer model.
"""
const qCalc = parsed_args["q"]

##################################################
# Options to be removed in GitHub
##################################################

"""
    gAni_default

Anisotropy parameter gAni for the Isochrone model.
"""
const gAni_default = parsed_args["gAni"]

"""
    logJrcutoff

log10(Jr) lower cut-off used in the computation of dRc/dt
"""
const logJrcutoff = parsed_args["logJrcutoff"]

"""
    logLcutoff

log10(L) lower cut-off used in the computation of dRc/dt
"""
const logLcutoff = parsed_args["logLcutoff"]

##################################################
##################################################

include("Constants.jl")
include("Mean.jl")
include("Bath.jl")
include("LocalDeflection.jl")
include("OrbitAverage.jl")
include("ActionSpace.jl")
include("Flux.jl")
