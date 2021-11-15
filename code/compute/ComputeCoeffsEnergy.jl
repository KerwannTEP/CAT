##################################################
# Compute the NR local diffusion coefficients in (E,L)
##################################################

include("../sources/julia/Main.jl")

r = 1.0 # Radius
E = -0.5 # Energy
L = 0.1 # Angular momentum
m = 0.1 # Mass of field star

println("--------------------")

println("Plummer potential")

println("q = ", qCalc)
println("m = ", m)
println("r = ", r)
println("E = ", E)
println("L = ", L)

println("--------------------")

println("Lc(E) = ", Lc(_tE(E)))
println("Coulomb log = ", logCoulomb)


println("--------------------")

@time dE, dE2, dL, dEdL, dL2 = localOrbitChange(r,E,L,m,100,m)


println("D_E  =  ", dE)
println("D_EE =  ", dE2)
println("D_L  =  ", dL)
println("D_LL =  ", dEdL)
println("D_EL =  ", dL2)

println("--------------------")
