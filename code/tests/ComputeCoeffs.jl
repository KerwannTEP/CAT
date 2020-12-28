include("../sources/julia/Main.jl")

E = 0.5
L = 0.1
q = 0.0
m = 0.1

@time dE, dE2, dL, dL2, dEdL = averageDiffCoeffs(E,L,q,m)
println("D_E  =  ",dE)
println("D_EE =  ",dE2)
println("D_L  =  ",dL)
println("D_LL =  ",dL2)
println("D_EL =  ",dEdL)