##################################################
# Compute a map of the NR diffusion coefficients 
##################################################

q_aniso  = 2.0
m_field  = 0.1
PARALLEL = "no"

##################################################

using HDF5

########################################
include("../sources/julia/Main.jl") # Loading the main code
########################################

EminMeasure, EmaxMeasure = 0.01,0.999
LminMeasure, LmaxMeasure = 0.01,2.5#Lc(EMeasure)-0.0001
nbEMeasure = 100
nbLMeasure = 100
#tabLMeasure = exp.(range(log(LminMeasure),length=nbLMeasure,log(LmaxMeasure)))
tabEMeasure = collect(range(EminMeasure,length=nbEMeasure,EmaxMeasure))
tabLMeasure = collect(range(LminMeasure,length=nbLMeasure,LmaxMeasure))

nbELGrid = nbEMeasure*nbLMeasure

const tabELGrid  = zeros(Float64,2,nbELGrid) # Location (E,L) of the grid points where the diffusion coefficients are computed

const tabDNRE  = zeros(Float64,nbELGrid) # Values of the DRR_E  coefficients on the (a,j)-grid
const tabDNREE = zeros(Float64,nbELGrid) # Values of the DRR_EE coefficients on the (a,j)-grid
const tabDNRL  = zeros(Float64,nbELGrid) # Values of the DRR_L  coefficients on the (a,j)-grid
const tabDNRLL = zeros(Float64,nbELGrid) # Values of the DRR_LL coefficients on the (a,j)-grid
const tabDNREL = zeros(Float64,nbELGrid) # Values of the DRR_EL coefficients on the (a,j)-grid

########################################

function tabELGrid!()
    index = 1
    for iE=1:nbEMeasure
    EMeasure = tabEMeasure[iE]
        for iL=1:nbLMeasure
            LMeasure = tabLMeasure[iL]
            tabELGrid[1,index], tabELGrid[2,index] = EMeasure, LMeasure
            index += 1
        end
    end
end

function tabDNR!()
    if (PARALLEL == "yes") 
        tid = 0
        PlummerTable_parallel = Array{IntTable,1}(undef,Threads.nthreads())
        for tid=1:Threads.nthreads()
            PlummerTable_parallel[tid] = IntTable_create!()
        end
        Threads.@threads for iGrid=1:nbELGrid
            EMeasure, LMeasure = tabELGrid[1,iGrid], tabELGrid[2,iGrid] 
            tid = Threads.threadid()

            if (EMeasure <= Ec(LMeasure))

                averageDiffCoeffs!(EMeasure,LMeasure,q_aniso,m_field,PlummerTable_parallel[tid])
                DE = PlummerTable_parallel[tid].dE[]
                DEE = PlummerTable_parallel[tid].dE2[]
                DL = PlummerTable_parallel[tid].dL[]
                DLL = PlummerTable_parallel[tid].dL2[]
                DEL = PlummerTable_parallel[tid].dEdL[]            
            else
                DE, DEE, DL, DLL, DEL = 0.0, 0.0, 0.0, 0.0, 0.0
            end
            tabDNRE[iGrid] = DE 
            tabDNREE[iGrid] = DEE 
            tabDNRL[iGrid] = DL 
            tabDNRLL[iGrid] = DLL 
            tabDNREL[iGrid] = DEL

        end
    else # Computation is not made in parallel
        for iGrid=1:nbELGrid
            EMeasure, LMeasure = tabELGrid[1,iGrid], tabELGrid[2,iGrid] 

            if (EMeasure <= Ec(LMeasure))

                averageDiffCoeffs!(EMeasure,LMeasure,q_aniso,m_field,PlummerTable_serial)
                DE = PlummerTable_serial.dE[]
                DEE = PlummerTable_serial.dE2[]
                DL = PlummerTable_serial.dL[]
                DLL = PlummerTable_serial.dL2[]
                DEL = PlummerTable_serial.dEdL[]            
            else
                DE, DEE, DL, DLL, DEL = 0.0, 0.0, 0.0, 0.0, 0.0
            end
            tabDNRE[iGrid] = DE 
            tabDNREE[iGrid] = DEE 
            tabDNRL[iGrid] = DL 
            tabDNRLL[iGrid] = DLL 
            tabDNREL[iGrid] = DEL

        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients_Orbital_Map_q_2.0.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabL",tabLMeasure) 
    write(file,"tabE",tabEMeasure) 
    write(file,"tabEL",tabELGrid) 
    write(file,"tabDNRE",tabDNRE) 
    write(file,"tabDNREE",tabDNREE) 
    write(file,"tabDNRL",tabDNRL) 
    write(file,"tabDNRLL",tabDNRLL) 
    write(file,"tabDNREL",tabDNREL) 
    close(file) # Closing the file
end
########################################

@time tabELGrid!()
@time tabDNR!()

########################################
writedump!(namefile) # Dumping the computed table