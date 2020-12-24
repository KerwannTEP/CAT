##################################################
# Useful general functions
##################################################

##################################################
# Compute Ec(L) the binding energy per unit mass of a circular orbit
##################################################

# Potential
function psi(r::Float64) 
    return 1.0/sqrt(1.0+r^2)
end

# Effective potential
function psiEff(r::Float64, L::Float64)
    return psi(r) - L^2/(2*r^2)
end

function dpsidr(r::Float64) 
    return -r/sqrt((1.0+r^2)^3)
end

function dpsiEffdr(r::Float64, L::Float64)
    return dpsidr(r) + L^2/(r^3)
end

function d2psidr2(r::Float64) 
    return -(sqrt((1.0+r^2)^3) - 3*r^2*sqrt(1.0+r^2))/(1.0+r^2)^3
end

function d2psiEffdr2(r::Float64, L::Float64)
    return d2psidr2(r) - 3* L^2/(r^4) 
end

# Finding the maximum of the effective potential with a Newton method applied on its derivative
function newton(L::Float64, eps::Float64=4.0*eps(Float64))
    r = L^(2/3)
    while(dpsidr(r+eps,L) > 0)
        r = r - dpsiEffdr(r,L)/d2psiEffdr2(r,L)
    end
    return r+eps/2
end

# Binding energy (per unit mass) for the circular orbit given an angular momentum L (per unit mass)
function Ec(L::Float64)
    if (L == 0.0)
        return 1.0
    end
    return psiEff(newton(L),L)
end

##################################################
# Convert (vr,vt) to (E,L) at a given radius r
##################################################

function bindingEnergy(vr::Float64, vt::Float64, r::Float64)
    return psi(r) - vr^2/(2.0) - vt^2/(2.0)
end

function angularMomentum(vr::Float64, vt::Float64, r::Float64)
    return r*vt
end

##################################################
# Convert (E,L) to (vr,vt) at a given radius r
# Convention vr >= 0
##################################################

function radialVelocity(E::Float64, L::Float64, r::Float64)
    return sqrt(2*(psiEff(r,L) - E))
end

function tangentVelocity(E::Float64, L::Float64, r::Float64)
    return L/r
end