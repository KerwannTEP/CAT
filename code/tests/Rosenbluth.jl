using HCubature

r = 1.0
E = 0.1
L = 0.1

function psi(r::Float64) 
    return 1.0/sqrt(1.0+r^2)
end

function psiEff(r::Float64, L::Float64)
    return psi(r) - L^2/(2*r^2)
end

function radialVel(r::Float64, E::Float64, L::Float64)
    return sqrt(2*abs((psiEff(r,L) - E)))
end

function tangentVel(r::Float64,E::Float64, L::Float64)
    return L/r
end

vr = radialVel(r,E,L)
vt = tangentVel(r,E,L)
v = sqrt(vr^2+vt^2)

function bindEnergy(r::Float64, vr::Float64, vt::Float64)
    return psi(r) - (vr^2)/2 - (vt^2)/2
end

function angMoment(r::Float64, vr::Float64, vt::Float64)
    return r*vt
end

function normVelIso(r::Float64, E::Float64)
    return psi(r) - (v^2)/2
end

function _vmax(r::Float64, vr::Float64, vt::Float64)
    return vr+vt+sqrt(2*(vr*vt+psi(r)))
end

vmax = _vmax(r,vr,vt)
MAXEVALS = 5000

function DF(E::Float64)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else # If E and L are positive
        return 3.0/(7.0*pi^3) * (2.0*E)^(7/2)
    end
end

function _hIso(r::Float64, v::Float64)
    int1 = (Ea->normVelIso(r,Ea[1])*DF(Ea[1]))
    int2 = (Ea->DF(Ea[1]))
    I1 = (1/v) * hcubature(int1,[E],[psi(r)],maxevals=MAXEVALS)[1]
    I2 = hcubature(int2,[0],[E],maxevals=MAXEVALS)[1]
    return 4*pi*(I1+I2)
end

println(_hIso(r,v))
