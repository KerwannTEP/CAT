include("../sources/julia/Main.jl")

r = 1.0
E = -0.5
L = 0.1
m = 0.1

println("--------------------")

println("Plummer potential")

println("q = ", qCalc)
println("E = ", E)
println("L = ", L)

println("--------------------")

println("Lc(E) = ", Lc(_tE(E)))
println("Coulomb log = ", logCoulomb)


println("--------------------")

@time dE, dE2, dL, dEdL, dL2 = localOrbitChange(r,E,L,m,qCalc,100,m)


println("D_E  =  ", dE)
println("D_EE =  ", dE2)
println("D_L  =  ", dL)
println("D_LL =  ", dEdL)
println("D_EL =  ", dL2)

println("--------------------")
