##################################################
# Compute a map of the NR diffusion coefficients
# In log-log action space
##################################################
using HDF5
########################################
include("../sources/julia/Main.jl") # Loading the main code
########################################

m_field  = _M/nbGlobularCluster

JrminMeasure, JrmaxMeasure = 0.00001*_L0,100.0*_L0
LminMeasure, LmaxMeasure = 0.00001*_L0,100.0*_L0
nbJrMeasure = 150
nbLMeasure = 150

tabJrMeasure = exp.(range(log(JrminMeasure),length=nbJrMeasure,log(JrmaxMeasure)))
tabLMeasure = exp.(range(log(LminMeasure),length=nbLMeasure,log(LmaxMeasure)))

nbJrLGrid = nbJrMeasure*nbLMeasure

const tabLJrGrid  = zeros(Float64,2,nbJrLGrid) # Location (E,L) of the grid points where the diffusion coefficients are computed

const tabRedDF  = zeros(Float64,nbJrLGrid) # Values of the DRR_E  coefficients on the (a,j)-grid

const tabDF  = zeros(Float64,nbJrLGrid) # Values of the DRR_E  coefficients on the (a,j)-grid
const tabDNRJr  = zeros(Float64,nbJrLGrid) # Values of the DRR_E  coefficients on the (a,j)-grid
const tabDNRJrTimeL  = zeros(Float64,nbJrLGrid) # Values of the DRR_L*L  coefficients on the (a,j)-grid
const tabDNRJrJr = zeros(Float64,nbJrLGrid) # Values of the DRR_EE coefficients on the (a,j)-grid
const tabDNRL  = zeros(Float64,nbJrLGrid) # Values of the DRR_L  coefficients on the (a,j)-grid
const tabDNRLTimeL  = zeros(Float64,nbJrLGrid) # Values of the DRR_L*L  coefficients on the (a,j)-grid
const tabDNRLL = zeros(Float64,nbJrLGrid) # Values of the DRR_LL coefficients on the (a,j)-grid
const tabDNRJrL = zeros(Float64,nbJrLGrid) # Values of the DRR_EL coefficients on the (a,j)-grid

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
            tabDNRJrTimeL[iGrid] = avrDJr*LMeasure
            tabDNRJrJr[iGrid] = avrDJrJr
            tabDNRL[iGrid] = avrDL
            tabDNRLTimeL[iGrid] = avrDL*LMeasure
            tabDNRLL[iGrid] = avrDLL
            tabDNRJrL[iGrid] = avrDJrL
        #    tabdEdJr[iGrid] = 1/dJrdE
    #        tabdLdJr[iGrid] = 1/dJrdL
        #    tabtE[iGrid] = tE
        #    tabtL[iGrid] = tL

        #    EMeasure = _E_from_Jr(JrMeasure,LMeasure)
        #    tabDF[iGrid] = _F(EMeasure,LMeasure)
        #    tabRedDF[iGrid] = 2*LMeasure*tabDF[iGrid]

        end
    else # Computation is not made in parallel
        for iGrid=1:nbJrLGrid
            LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]

            avrDJr, avrDL, avrDJrJr, avrDJrL, avrDLL = avr_action_coefficients(JrMeasure,LMeasure,m_field)

            tabDNRJr[iGrid] = avrDJr
            tabDNRJrTimeL[iGrid] = avrDJr*LMeasure
            tabDNRJrJr[iGrid] = avrDJrJr
            tabDNRL[iGrid] = avrDL
            tabDNRLTimeL[iGrid] = avrDL*LMeasure
            tabDNRLL[iGrid] = avrDLL
            tabDNRJrL[iGrid] = avrDJrL
        #    tabdEdJr[iGrid] = 1/dJrdE
        #    tabdLdJr[iGrid] = 1/dJrdL
    #        tabtE[iGrid] = tE
    #        tabtL[iGrid] = tL

    #        EMeasure = _E_from_Jr(JrMeasure,LMeasure)
        #    tabDF[iGrid] = _F(EMeasure,LMeasure)
        #    tabRedDF[iGrid] = 2*LMeasure*tabDF[iGrid]

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
    write(file,"tabDNRJrTimeL",tabDNRJrTimeL)
    write(file,"tabDNRLTimeL",tabDNRLTimeL)
    write(file,"tabDNRLL",tabDNRLL)
    write(file,"tabDNRJrL",tabDNRJrL)
#    write(file,"tabdEdJr",tabdEdJr)
#    write(file,"tabdLdJr",tabdLdJr)
#    write(file,"tabtE",tabtE)
#    write(file,"tabtL",tabtL)
#    write(file,"tabDF",tabDF)
#    write(file,"tabRedDF",tabRedDF)
    write(file,"nbJrMeasure",nbJrMeasure)
    write(file,"nbLMeasure",nbLMeasure)


    close(file) # Closing the file
end
########################################

@time tabLJrGrid!()
@time tabDNR!()

########################################
writedump!(namefile) # Dumping the computed table
