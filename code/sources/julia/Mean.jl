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
function maximizerPsiEff(L::Float64, eps::Float64=4.0*eps(Float64))
    r = L^(2/3)
    while(dpsiEffdr(r+eps,L) > 0)
        r = r - dpsiEffdr(r,L)/d2psiEffdr2(r,L)
    end
    return r+eps/2
end

# Binding energy (per unit mass) for the circular orbit given an angular momentum L (per unit mass)
function Ec(L::Float64)
    if (L == 0.0)
        return 1.0
    else
        return psiEff(maximizerPsiEff(L),L)
    end
end

##################################################
# Compute Lc(E) the angular momentum per unit mass of a circular orbit
##################################################
# Issue

function _z(r::Float64, E::Float64)
    return 2*r^2*(psi(r)-E)
end

function _dzdr(r::Float64, E::Float64)
    return 4*r*(psi(r)-E) + 2*r^2*dpsidr(r)
end

function _d2zdr2(r::Float64, E::Float64)
    return 4*(psi(r)-E) + 8*r*dpsidr(r) + 2*r^2*d2psidr2(r)
end

function maximizerZ(E::Float64, eps::Float64=4.0*eps(Float64))
    r = sqrt(E^(-2)-1)
    while(_dzdr(r-eps,E) < 0)
        r = r - _dzdr(r,E)/_d2zdr2(r,E)
    end
    return r-eps/2
end

function Lc(E::Float64)
    if (E == 1.0)
        return 0.0
    else
        return sqrt(_z(maximizerZ(E),E))
    end
end

##################################################
# Convert (vr,vt) to (E,L) at a given radius r
##################################################

function bindingEnergy(r::Float64, vr::Float64, vt::Float64)
    return psi(r) - vr^2/(2.0) - vt^2/(2.0)
end

function angularMomentum(r::Float64, vr::Float64, vt::Float64)
    return r*vt
end

##################################################
# Convert (E,L) to (vr,vt) at a given radius r
# Convention vr >= 0
##################################################

function radialVelocity(r::Float64,E::Float64, L::Float64)
    return sqrt(2*abs((psiEff(r,L) - E)))
end

function tangentVelocity(r::Float64,E::Float64, L::Float64)
    return L/r
end

##################################################
# Solve the degree-3 polynomial equation aX^3 + bX^2 + cX + d = 0
# Gives the real roots
##################################################

function solveRealPoly3(a::Float64, b::Float64, c::Float64, d::Float64)
    p = (3*a*c - b^2)/(3*a^2)
    q = (2*b^3 - 9*a*b*c + 27*a^2*d)/(27*a^3)
    disc = -(4*p^3 + 27*q^2)
    if (disc < 0.0)
        rt1 = cbrt(-q/2 + sqrt((q^2)/4 + (p^3)/27)) + cbrt(-q/2 - sqrt((q^2)/4 + (p^3)/27))
        rt1 -= b/(3*a)
        nbRoot = 1
        return nbRoot, rt1
    else
        rt1 = 2*sqrt(-p/3)*cos((1/3)*acos((3*q)/(2*p) * sqrt(-3/p)))
        rt2 = 2*sqrt(-p/3)*cos((1/3)*acos((3*q)/(2*p) * sqrt(-3/p)) - 2*pi/3)
        rt3 = 2*sqrt(-p/3)*cos((1/3)*acos((3*q)/(2*p) * sqrt(-3/p)) - 4*pi/3) 
        rt1 -= b/(3*a)
        rt2 -= b/(3*a)
        rt3 -= b/(3*a)
        nbRoot = 3
        return nbRoot, rt1, rt2, rt3
    end
end