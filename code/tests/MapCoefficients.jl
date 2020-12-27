##################################################
# Compute a map of the NR diffusion coefficients at some fixed L
# For E between Emin and Ec(L)
##################################################

LMeasure = 0.1
q_aniso  = 0.0
m_field  = 0.1
PARALLEL = "no"

##################################################

using HDF5

########################################
include("../sources/Main.jl") # Loading the main code
########################################

EminMeasure, EmaxMeasure = 0.01, Ec(L)-0.0001
nbEMeasure = 30
tabEMeasure = exp.(range(log(EminMeasure),length=nbEMeasure,log(EmaxMeasure)))

const tabDNRE  = zeros(Float64,nbjMeasure) # Values of the DRR_E  coefficients on the (a,j)-grid
const tabDNREE = zeros(Float64,nbjMeasure) # Values of the DRR_EE coefficients on the (a,j)-grid
const tabDNRL  = zeros(Float64,nbjMeasure) # Values of the DRR_L  coefficients on the (a,j)-grid
const tabDNRLL = zeros(Float64,nbjMeasure) # Values of the DRR_LL coefficients on the (a,j)-grid
const tabDNREL = zeros(Float64,nbjMeasure) # Values of the DRR_EL coefficients on the (a,j)-grid

function tabDNR!()
    if (PARALLEL == "yes") 
        Threads.@threads for iE=1:nbEMeasure 
            EMeasure = tabEMeasure[iE] 
            DE, DEE, DL, DLL, DEL = averageDiffCoeffs(EMeasure,LMeasure,q_aniso,m_field)

            tabDNRE[iE] = DE 
            tabDNREE[iE] = DEE 
            tabDNRL[iE] = DL 
            tabDNRLL[iE] = DLL 
            tabDNREL[iE] = DEL

        end
    else # Computation is not made in parallel
        for iE=1:nbEMeasure 
            EMeasure = tabEMeasure[iE] 
            DE, DEE, DL, DLL, DEL = averageDiffCoeffs(EMeasure,LMeasure,q_aniso,m_field)

            tabDNRE[iE] = DE 
            tabDNREE[iE] = DEE 
            tabDNRL[iE] = DL 
            tabDNRLL[iE] = DLL 
            tabDNREL[iE] = DEL

        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients_Fixed_L.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabE",tabEMeasure) 
    write(file,"LMeasure",LMeasure) 
    write(file,"tabDNRE",tabDNRE) 
    write(file,"tabDNREE",tabDNREE) 
    write(file,"tabDNRL",tabDNRL) 
    write(file,"tabDNRLL",tabDNRLL) 
    write(file,"tabDNREL",tabDNREL) 
    close(file) # Closing the file
end
########################################

@time tabDNR!()

########################################
writedump!(namefile) # Dumping the computed table