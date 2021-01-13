include("../sources/julia/Main.jl")

E = 0.7
L = 0.2
q = 0.0
m = 0.1

println("Anisotropic")
@time averageDiffCoeffs!(E,L,q,m)
aE = PlummerTable_serial.dE[]
aE2 = PlummerTable_serial.dE2[]
aL = PlummerTable_serial.dL[]
aL2 = PlummerTable_serial.dL2[]
aEL = PlummerTable_serial.dEdL[]
println("D_E  =  ",aE)
println("D_EE =  ",aE2)
println("D_L  =  ",aL)
println("D_LL =  ",aL2)
println("D_EL =  ",aEL)

println("--------------------")
println("Isotropic")

@time averageDiffCoeffsIso!(E,L,m)
iE = PlummerTable_serial.dE[]
iE2 = PlummerTable_serial.dE2[]
iL = PlummerTable_serial.dL[]
iL2 = PlummerTable_serial.dL2[]
iEL = PlummerTable_serial.dEdL[]
println("D_E  =  ",iE)
println("D_EE =  ",iE2)
println("D_L  =  ",iL)
println("D_LL =  ",iL2)
println("D_EL =  ",iEL)

println("--------------------")
println("Relative Error Aniso/Isotropic")

println("D_E  :  ",trunc(Int,abs(iE-aE)/abs(iE)*1000),"‰")
println("D_EE :  ",trunc(Int,abs(iE2-aE2)/abs(iE2)*1000),"‰")
println("D_L  :  ",trunc(Int,abs(iL-aL)/abs(iL)*1000),"‰")
println("D_LL :  ",trunc(Int,abs(iL2-aL2)/abs(iL2)*1000),"‰")
println("D_EL :  ",trunc(Int,abs(iEL-aEL)/abs(iEL)*1000),"‰")
