##################################################
# Orbit-averagering
##################################################

using QuadGK

##################################################
# Determine the nature of the orbit
# Use positive binding energy
##################################################

function _alpha(E::Float64, L::Float64)
    return 4.0*E^2
end

function _beta(E::Float64, L::Float64)
    return 4.0*(E^2-1.0+E*L^2)
end

function _gamma(E::Float64, L::Float64)
    return 4.0*E*L^2+L^4
end

function _delta(E::Float64, L::Float64)
    return L^4
end

function _p(E::Float64, L::Float64)
    a = _alpha(E,L)
    b = _beta(E,L)
    c = _gamma(E,L)
    return (3.0*a*c-b^2)/(3.0*a^2)
end

function _q(E::Float64, L::Float64)
    a = _alpha(E,L)
    b = _beta(E,L)
    c = _gamma(E,L)
    d = _delta(E,L)
    return (2.0*b^3-9.0*a*b*c+27.0*a^2*d)/(27.0*a^3)
end

function discriminantY(p::Float64, q::Float64)
    return -(4.0*p^3+27*q^2)
end

# Returns a tuple (number of roots, root1, ...) of the roots of the X-polynomial
# Returns (0,"Unbounded orbit") if E <= 0
function rootX(E::Float64, L::Float64)
    if (E <= 0.0)
        return (0,"Unbounded orbit")
    end
    a = _alpha(E,L)
    b = _beta(E,L)
    p = _p(E,L)
    q = _q(E,L)
    disc = discriminantY(p,q)
    if (disc < 0.0)
        rt = cbrt(-q/2+sqrt(q^2/4+p^3/27)) + cbrt(-q/2-sqrt(q^2/4+p^3/27))
        rt += b/(3.0*a)
        return (1,rt)
    end
    rt1 = 2*sqrt(-p/3.0)*cos(1.0/3 * acos(3.0*q/(2.0*p) * sqrt(-3.0/p)))
    rt1 -= b/(3.0*a)
    rt2 = 2*sqrt(-p/3.0)*cos(1.0/3 * acos(3.0*q/(2.0*p) * sqrt(-3.0/p))-2.0*pi/3)
    rt2 -= b/(3.0*a)
    rt3 = 2*sqrt(-p/3.0)*cos(1.0/3 * acos(3.0*q/(2.0*p) * sqrt(-3.0/p))-2.0*pi*2.0/3)
    rt3 -= b/(3.0*a)
    return (3, [rt1, rt2, rt3])
end

# Returns "Unbounded orbit" if E<=0 and (rmin,rmax) if the parameters are allowed
# Check beforehands that the parameters are allowed, i.e!
# -> E,L are positive
# -> E < Ec(L)
function radiusBounds(E::Float64, L::Float64)
    nbRoots, rts = rootX(E,L)
    if (nbRoots == 0)
        return rts
    end
    disc = discriminantY(E,L)

    # Check if the parameters (E,L) can describe an orbit
    @assert nbRoots == 3 "radiusBounds: Forbidden parameters (E,L)"
    
    rts = sort(rts)  # Sort the roots
    rmin, rmax = rts[2], rts[3] # Lowest root is regative which the two other are positive
    rmin = sqrt(rmin)
    rmax = sqrt(rmax)    
    return rmin, rmax
end

##################################################
# Compute the period of an orbit
##################################################

function periodOrbit(E::Float64, L::Float64)
    rmin, rmax = radiusBounds(E,L)
    period = quadgk(r->sqrt(2.0*(psiEff(r,L)-E)),rmin,rmax)[1]
    return period
end

##################################################
# Compute the orbit-averaged NR diffusion coefficients
# interpolate them beforehand
##################################################