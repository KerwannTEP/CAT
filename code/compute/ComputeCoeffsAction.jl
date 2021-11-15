##################################################
# Compute the NR orbit-averaged diffusion coefficients in (Jr,L)
##################################################

include("../sources/julia/Main.jl")

Jr = 0.2 # Radial action
L = 0.3 # Angular momentum
m = 0.1 # Mass of field star

println("--------------------")

println("Plummer potential")

println("q  = ", qCalc)
println("m  = ", m)
println("Jr = ", Jr)
println("L  = ", L)

println("--------------------")

println("Coulomb log = ", logCoulomb)


println("--------------------")

@time dJr, dL, dJrJr, dJrL, dLL = avr_action_coefficients(Jr,L,m)


println("D_Jr  =  ", dJr)
println("D_JrJr =  ", dJrJr)
println("D_L  =  ", dL)
println("D_LL =  ", dLL)
println("D_JrL =  ", dJrL)

println("--------------------")
