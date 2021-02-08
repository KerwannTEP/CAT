include("../sources/julia/Main.jl")

r = 0.1
E = 0.1
L = 0.1
q = 0.0
m = 0.1
vr = radialVelocity(r,E,L)
vt = tangentVelocity(r,E,L)
v = sqrt(vr^2+vt^2)

println("Anisotropic")
@time localVelChange!(r,vr,vt,q,m)
avp = PlummerTable_serial.dvPar[]
avp2= PlummerTable_serial.dvPar2[]
avt2 = PlummerTable_serial.dvTan2[]
println("dvPar  = ",avp)
println("dvPar2 = ",avp2)
println("dvTan2   ",avt2)

println("--------------------")
println("Isotropic")

@time localVelChangeIso!(r,v,m)
ivp = PlummerTable_serial.dvPar[]
ivp2= PlummerTable_serial.dvPar2[]
ivt2 = PlummerTable_serial.dvTan2[]
println("dvPar  = ",ivp)
println("dvPar2 = ",ivp2)
println("dvTan2   ",ivt2)

println("--------------------")
println("Relative Error Aniso/Isotropic")

#println("dvPar  :  ",trunc(Int,abs(ivp-avp)/abs(ivp)*1000),"‰")
#println("dvPar2 :  ",trunc(Int,abs(ivp2-avp2)/abs(ivp2)*1000),"‰")
#println("dvTan2 :  ",trunc(Int,abs(ivt2-avt2)/abs(ivt2)*1000),"‰")

println("dvPar  :  ",abs(ivp-avp)/abs(ivp))
println("dvPar2 :  ",abs(ivp2-avp2)/abs(ivp2))
println("dvTan2 :  ",abs(ivt2-avt2)/abs(ivt2))
