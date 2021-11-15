##################################################
# Compute a map of the dF/dt in action space
# Linear sampling
##################################################
using HDF5
########################################
include("../sources/julia/Main.jl") # Loading the main code
include("../sources/julia/test/OldFunctions.jl")
########################################

m_field  = _M/nbGlobularCluster # Mass of field star
epsRef = _L0*10^(-5)

########################################
# Action space parameter
########################################

JrminMeasure, JrmaxMeasure = _L0*0.0001, _L0*0.5 # Jr range
LminMeasure, LmaxMeasure = _L0*0.0001, _L0*1.0 # L range
nbJrMeasure = 100 # Number of Jr sampling points
nbLMeasure = 500 # Number of L sampling points

########################################

tabJrMeasure = collect(range(JrminMeasure,length=nbJrMeasure,JrmaxMeasure))
tabLMeasure = collect(range(LminMeasure,length=nbLMeasure,LmaxMeasure))
nbJrLGrid = nbJrMeasure*nbLMeasure

const tabLJrGrid  = zeros(Float64,2,nbJrLGrid) # Location (L,Jr) of the grid points where the diffusion coefficients are computed
const tabdFdt  = zeros(Float64,nbJrLGrid) # Values of dF/dt on the (L,Jr)-grid


########################################
# Functions to fill the arrays of coefficients
########################################

function tabLJrGrid!()
    index = 1
    for iJr=1:nbJrMeasure
    JrMeasure = tabJrMeasure[iJr]
        for iL=1:nbLMeasure
            LMeasure = tabLMeasure[iL]
            tabLJrGrid[1,index], tabLJrGrid[2,index] = LMeasure, JrMeasure
            index += 1
        end
    end
end

function tabdFdt!()
    if (PARALLEL == "yes")
        println("parallel")
        Threads.@threads for iGrid=1:nbJrLGrid
            LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]
            eps = min(epsRef,JrMeasure*10^(-3),LMeasure*10^(-3))

            dfdt = -divflux_NR_num_expand(JrMeasure,LMeasure,m_field,eps)
            tabdFdt[iGrid] = dfdt
        end
    else # Computation is not made in parallel
        for iGrid=1:nbJrLGrid
            LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]
            eps = min(epsRef,JrMeasure*10^(-3),LMeasure*10^(-3))

            dfdt = -divflux_NR_num_expand(JrMeasure,LMeasure,m_field,eps)
            tabdFdt[iGrid] = dfdt
        end
    end
end

########################################
namefile = "../data/Dump_dFdt_Map_q_"*string(qCalc)*".hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"q",qCalc)
    write(file,"tabJr",tabJrMeasure)
    write(file,"tabL",tabLMeasure)
    write(file,"tabLJr",tabLJrGrid)
    write(file,"tabdFdt",tabdFdt)

    write(file,"nbJrMeasure",nbJrMeasure)
    write(file,"nbLMeasure",nbLMeasure)
    write(file,"Mtot",_M)
    write(file,"Npart",nbGlobularCluster)

    close(file) # Closing the file
end
########################################

@time tabLJrGrid!()
@time tabdFdt!()

########################################
writedump!(namefile) # Dumping the computed table
