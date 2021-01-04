##################################################
# Compute a map of the NR diffusion coefficients at some fixed E
# For L between Emin and Emax
##################################################

EMeasure = 0.2
q_aniso  = -10.1
m_field  = 0.1
PARALLEL = "no"

##################################################

using HDF5

########################################
include("../sources/julia/Main.jl") # Loading the main code
########################################

LminMeasure, LmaxMeasure = 0.01,1.5
nbLMeasure = 20
#tabLMeasure = exp.(range(log(LminMeasure),length=nbLMeasure,log(LmaxMeasure)))
tabLMeasure = collect(range(LminMeasure,length=nbLMeasure,LmaxMeasure))

const tabDNRE  = zeros(Float64,nbLMeasure) # Values of the DRR_E  coefficients on the (a,j)-grid
const tabDNREE = zeros(Float64,nbLMeasure) # Values of the DRR_EE coefficients on the (a,j)-grid
const tabDNRL  = zeros(Float64,nbLMeasure) # Values of the DRR_L  coefficients on the (a,j)-grid
const tabDNRLL = zeros(Float64,nbLMeasure) # Values of the DRR_LL coefficients on the (a,j)-grid
const tabDNREL = zeros(Float64,nbLMeasure) # Values of the DRR_EL coefficients on the (a,j)-grid

function tabDNR!()
    if (PARALLEL == "yes") 
        tid = 0
        PlummerTable_parallel = Array{IntTable,1}(undef,Threads.nthreads())
        for tid=1:Threads.nthreads()
            PlummerTable_parallel[tid] = IntTable_create!()
        end
        Threads.@threads for iL=1:nbLMeasure 
            LMeasure = tabLMeasure[iL] 
            tid = Threads.threadid()
#println(tid)
#println(PlummerTable_parallel[tid])
            
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
                averageDiffCoeffs!(EMeasure,LMeasure,q_aniso,m_field,PlummerTable_serial)
                DE = PlummerTable_serial.dE[]
                DEE = PlummerTable_serial.dE2[]
                DL = PlummerTable_serial.dL[]
                DLL = PlummerTable_serial.dL2[]
                DEL = PlummerTable_serial.dEdL[]
            else
                DE, DEE, DL, DLL, DEL = 0.0, 0.0, 0.0, 0.0, 0.0
            end
            tabDNRE[iL] = DE 
            tabDNREE[iL] = DEE 
            tabDNRL[iL] = DL 
            tabDNRLL[iL] = DLL 
            tabDNREL[iL] = DEL

            #GC.gc()
            
        end
    end
end

########################################
namefile = "../data/Dump_Diffusion_Coefficients_Fixed_E_0.2.hf5"
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