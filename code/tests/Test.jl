include("../sources/Main.jl")

E = 0.4
L = 0.1

p = _p(E,L)
q = _q(E,L)
@time disc = discriminantY(p,q)
println("Discriminant =  ",disc)

@time rts = radiusBounds(E,L)
println("Radius bounds : ",rts)