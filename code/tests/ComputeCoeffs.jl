include("../sources/julia/Main.jl")

E = 0.5
L = 0.1
q = 0.0
m = 0.1

@time averageDiffCoeffs!(E,L,q,m)
println("D_E  =  ",PlummerTable_serial.dE[])
println("D_EE =  ",PlummerTable_serial.dE2[])
println("D_L  =  ",PlummerTable_serial.dL[])
println("D_LL =  ",PlummerTable_serial.dL2[])
println("D_EL =  ",PlummerTable_serial.dEdL[])