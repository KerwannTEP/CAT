##################################################
# Compute a map of the NR diffusion coefficients
# In log-log action space
##################################################
using HDF5
########################################
include("../sources/julia/Main.jl") # Loading the main code
########################################

m_field  = _M/nbGlobularCluster # Mass of field star

########################################
# Action space parameter
########################################

JrminMeasure, JrmaxMeasure = 0.00001*_L0,100.0*_L0 # Jr range
LminMeasure, LmaxMeasure = 0.00001*_L0,100.0*_L0 # L range
nbJrMeasure = 150 # Number of Jr sampling points
nbLMeasure = 150 # Number of L sampling points

########################################

tabJrMeasure = exp.(range(log(JrminMeasure),length=nbJrMeasure,log(JrmaxMeasure)))
tabLMeasure = exp.(range(log(LminMeasure),length=nbLMeasure,log(LmaxMeasure)))
nbJrLGrid = nbJrMeasure*nbLMeasure

const tabLJrGrid  = zeros(Float64,2,nbJrLGrid) # Location (L,Jr) of the grid points where the diffusion coefficients are computed
const tabDNRJr  = zeros(Float64,nbJrLGrid) # Values of the DNR_Jr coefficients on the (L,Jr)-grid
const tabDNRJrJr = zeros(Float64,nbJrLGrid) # Values of the DNR_JrJr coefficients on the (L,Jr)-grid
const tabDNRL  = zeros(Float64,nbJrLGrid) # Values of the DNR_L coefficients on the (L,Jr)-grid
const tabDNRLL = zeros(Float64,nbJrLGrid) # Values of the DNR_LL coefficients on the (L,Jr)-grid
const tabDNRJrL = zeros(Float64,nbJrLGrid) # Values of the DNR_JrL coefficients on the (L,Jr)-grid

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

function tabDNR!()
    if (PARALLEL == "yes")
        println("parallel")
        Threads.@threads for iGrid=1:nbJrLGrid
            LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]

            avrDJr, avrDL, avrDJrJr, avrDJrL, avrDLL = avr_action_coefficients(JrMeasure,LMeasure,m_field)

            tabDNRJr[iGrid] = avrDJr
            tabDNRJrJr[iGrid] = avrDJrJr
            tabDNRL[iGrid] = avrDL
            tabDNRLL[iGrid] = avrDLL
            tabDNRJrL[iGrid] = avrDJrL

        end
    else # Computation is not made in parallel
        for iGrid=1:nbJrLGrid
            LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]

            avrDJr, avrDL, avrDJrJr, avrDJrL, avrDLL = avr_action_coefficients(JrMeasure,LMeasure,m_field)

            tabDNRJr[iGrid] = avrDJr
            tabDNRJrJr[iGrid] = avrDJrJr
            tabDNRL[iGrid] = avrDL
            tabDNRLL[iGrid] = avrDLL
            tabDNRJrL[iGrid] = avrDJrL

        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients_Action_Orbital_Map_q_"*string(qCalc)*".hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabL",tabLMeasure)
    write(file,"q",qCalc)
    write(file,"tabJr",tabJrMeasure)
    write(file,"tabLJr",tabLJrGrid)
    write(file,"tabDNRJr",tabDNRJr)
    write(file,"tabDNRJrJr",tabDNRJrJr)
    write(file,"tabDNRL",tabDNRL)
    write(file,"tabDNRLL",tabDNRLL)
    write(file,"tabDNRJrL",tabDNRJrL)

    write(file,"nbJrMeasure",nbJrMeasure)
    write(file,"nbLMeasure",nbLMeasure)

    close(file) # Closing the file
end
########################################

@time tabLJrGrid!()
@time tabDNR!()

########################################
writedump!(namefile) # Dumping the computed table
