##################################################
# Compute a map of the dF/dt
##################################################


##################################################

using HDF5

########################################
include("../sources/julia/Main.jl") # Loading the main code
include("../sources/julia/test/OldFunctions.jl")



########################################

#Npart = 10^5
m_field  = _M/nbGlobularCluster
epsRef = _L0*10^(-5)

JrminMeasure, JrmaxMeasure = _L0*0.0001, _L0*0.5 #0.5
LminMeasure, LmaxMeasure = _L0*0.0001, _L0*1.0
nbJrMeasure = 100
nbLMeasure = 500

tabJrMeasure = collect(range(JrminMeasure,length=nbJrMeasure,JrmaxMeasure))
tabLMeasure = collect(range(LminMeasure,length=nbLMeasure,LmaxMeasure))

# tabJrMeasure = exp.(range(log(JrminMeasure),length=nbJrMeasure,log(JrmaxMeasure)))
# tabLMeasure = exp.(range(log(LminMeasure),length=nbLMeasure,log(LmaxMeasure)))


println(tabJrMeasure)
println(tabLMeasure)

nbJrLGrid = nbJrMeasure*nbLMeasure

const tabLJrGrid  = zeros(Float64,2,nbJrLGrid)

const tabdFdt  = zeros(Float64,nbJrLGrid)

const tabDF  = zeros(Float64,nbJrLGrid)
const tabdDFdE  = zeros(Float64,nbJrLGrid)
const tabdDFdL  = zeros(Float64,nbJrLGrid)
const tabdDFdEE  = zeros(Float64,nbJrLGrid)
const tabdDFdEL  = zeros(Float64,nbJrLGrid)
const tabdDFdLL  = zeros(Float64,nbJrLGrid)
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
            EMeasure = _E_from_Jr(JrMeasure,LMeasure)
            Ftot = _F(EMeasure,LMeasure)
            redDF = 2*LMeasure*Ftot

            dE, dL, dEE, dEL, dLL = _dF(EMeasure,LMeasure)
            dFdE = 2*LMeasure*dE
            dFdL = 2*Ftot+2*LMeasure*dL
            dFdEE = 2*LMeasure*dEE
            dFdEL = 2*dE+2*LMeasure*dEL
            dFdLL = 2*dL+2*dL+2*LMeasure*dLL

            tabdFdt[iGrid] = dfdt
            tabDF[iGrid] = redDF

            tabdDFdE[iGrid] = dFdE
            tabdDFdL[iGrid] = dFdL
            tabdDFdEE[iGrid] = dFdEE
            tabdDFdEL[iGrid] = dFdEL
            tabdDFdLL[iGrid] = dFdLL


        end
    else # Computation is not made in parallel
        for iGrid=1:nbJrLGrid
            LMeasure, JrMeasure = tabLJrGrid[1,iGrid], tabLJrGrid[2,iGrid]

            eps = min(epsRef,JrMeasure*10^(-3),LMeasure*10^(-3))

            dfdt = -divflux_NR_num_expand(JrMeasure,LMeasure,m_field,eps)
            EMeasure = _E_from_Jr(JrMeasure,LMeasure)
            Ftot = _F(EMeasure,LMeasure)
            redDF = 2*LMeasure*Ftot

            dE, dL, dEE, dEL, dLL = _dF(EMeasure,LMeasure)
            dFdE = 2*LMeasure*dE
            dFdL = 2*Ftot+2*LMeasure*dL
            dFdEE = 2*LMeasure*dEE
            dFdEL = 2*dE+2*LMeasure*dEL
            dFdLL = 2*dL+2*dL+2*LMeasure*dLL

            tabdFdt[iGrid] = dfdt
            tabDF[iGrid] = redDF

            tabdDFdE[iGrid] = dFdE
            tabdDFdL[iGrid] = dFdL
            tabdDFdEE[iGrid] = dFdEE
            tabdDFdEL[iGrid] = dFdEL
            tabdDFdLL[iGrid] = dFdLL

        end
    end
end

########################################
#namefile = "../data/Dump_dFdt_Map_q_"*string(qCalc)*".hf5"
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
    write(file,"tabDF",tabDF)
    write(file,"nbJrMeasure",nbJrMeasure)
    write(file,"nbLMeasure",nbLMeasure)
    write(file,"Mtot",_M)
    write(file,"Npart",nbGlobularCluster)

    write(file,"tabdDFdE",tabdDFdE)
    write(file,"tabdDFdL",tabdDFdL)
    write(file,"tabdDFdEE",tabdDFdEE)
    write(file,"tabdDFdEL",tabdDFdEL)
    write(file,"tabdDFdLL",tabdDFdLL)

    close(file) # Closing the file
end
########################################

@time tabLJrGrid!()
@time tabdFdt!()

########################################
writedump!(namefile) # Dumping the computed table
