##################################################
# Compute a map of the NR diffusion coefficients at some fixed E
# For L between Emin and Emax
##################################################

EMeasure = 0.1
q_aniso  = 0.0
m_field  = 0.1
PARALLEL = "no"

##################################################

using HDF5

########################################
include("../sources/Main.jl") # Loading the main code
########################################

LminMeasure, LmaxMeasure = 0.001,2.2
nbLMeasure = 200
tabLMeasure = exp.(range(log(LminMeasure),length=nbLMeasure,log(LmaxMeasure)))

const tabDNRE  = zeros(Float64,nbLMeasure) # Values of the DRR_E  coefficients on the (a,j)-grid
const tabDNREE = zeros(Float64,nbLMeasure) # Values of the DRR_EE coefficients on the (a,j)-grid
const tabDNRL  = zeros(Float64,nbLMeasure) # Values of the DRR_L  coefficients on the (a,j)-grid
const tabDNRLL = zeros(Float64,nbLMeasure) # Values of the DRR_LL coefficients on the (a,j)-grid
const tabDNREL = zeros(Float64,nbLMeasure) # Values of the DRR_EL coefficients on the (a,j)-grid

function tabDNR!()
    if (PARALLEL == "yes") 
        Threads.@threads for iL=1:nbLMeasure 
            LMeasure = tabLMeasure[iL] 
            if (EMeasure <= Ec(LMeasure))
                DE, DEE, DL, DLL, DEL = averageDiffCoeffs(EMeasure,LMeasure,q_aniso,m_field)
            else
                DE, DEE, DL, DLL, DEL = 0.0, 0.0, 0.0, 0.0, 0.0
            end
            tabDNRE[iL] = DE 
            tabDNREE[iL] = DEE 
            tabDNRL[iL] = DL 
            tabDNRLL[iL] = DLL 
            tabDNREL[iL] = DEL

        end
    else # Computation is not made in parallel
        for iL=1:nbLMeasure 
            LMeasure = tabLMeasure[iL] 
            if (EMeasure <= Ec(LMeasure))
                DE, DEE, DL, DLL, DEL = averageDiffCoeffs(EMeasure,LMeasure,q_aniso,m_field)
            else
                DE, DEE, DL, DLL, DEL = 0.0, 0.0, 0.0, 0.0, 0.0
            end
            tabDNRE[iL] = DE 
            tabDNREE[iL] = DEE 
            tabDNRL[iL] = DL 
            tabDNRLL[iL] = DLL 
            tabDNREL[iL] = DEL
            
        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients_Fixed_E.hf5"
########################################
# Function that writes to .hf5 files
########################################
function writedump!(namefile)
    file = h5open(namefile,"w") # Opening the file
    write(file,"tabL",tabLMeasure) 
    write(file,"EMeasure",EMeasure) 
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