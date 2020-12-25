##################################################
# Orbit-averagering
##################################################

using QuadGK
using Interpolations

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
    else
        a = _alpha(E,L)
        b = _beta(E,L)
        p = _p(E,L)
        q = _q(E,L)
        disc = discriminantY(p,q)
        if (disc < 0.0)
            rt = cbrt(-q/2+sqrt(q^2/4+p^3/27)) + cbrt(-q/2-sqrt(q^2/4+p^3/27))
            rt += b/(3.0*a)
            return (1,rt)
        else
            rt1 = 2*sqrt(-p/3.0)*cos(1.0/3 * acos(3.0*q/(2.0*p) * sqrt(-3.0/p)))
            rt1 -= b/(3.0*a)
            rt2 = 2*sqrt(-p/3.0)*cos(1.0/3 * acos(3.0*q/(2.0*p) * sqrt(-3.0/p))-2.0*pi/3)
            rt2 -= b/(3.0*a)
            rt3 = 2*sqrt(-p/3.0)*cos(1.0/3 * acos(3.0*q/(2.0*p) * sqrt(-3.0/p))-2.0*pi*2.0/3)
            rt3 -= b/(3.0*a)
            return (3, [rt1, rt2, rt3])
        end
    end
end

# Returns "Unbounded orbit" if E<=0 and (rmin,rmax) if the parameters are allowed
# Check beforehands that the parameters are allowed, i.e!
# -> E,L are positive
# -> E < Ec(L)
function radiusBounds(E::Float64, L::Float64)
    nbRoots, rts = rootX(E,L)
    if (nbRoots == 0)
        return rts
    else
        disc = discriminantY(E,L)
        # Check if the parameters (E,L) can describe an orbit
        @assert nbRoots == 3 "radiusBounds: Forbidden parameters (E,L)"
    
        rts = sort(rts)  # Sort the roots
        rmin, rmax = rts[2], rts[3] # Lowest root is regative which the two other are positive
        rmin = sqrt(rmin)
        rmax = sqrt(rmax)    
        return rmin, rmax
    end
end

##################################################
# Compute the period of an orbit
##################################################

function halfPeriodOrbit(E::Float64, L::Float64)
    rmin, rmax = radiusBounds(E,L)
    halfperiod = quadgk(r->sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
    return halfperiod
end

##################################################
# Control the error on the period due to uncertainty on the boundaries
##################################################

function halfPeriodUncertainty(E::Float64, L::Float64)
    rmin, rmax = radiusBounds(E,L)
    dt = sqrt(abs(psiEff(rmin,L)-E))/abs(dpsiEffdr(rmin,L)) + 
         sqrt(abs(psiEff(rmax,L)-E))/abs(dpsiEffdr(rmax,L))
    return dt
end

function relativeHalfPeriodUncertainty(E::Float64, L::Float64)
    t = halfPeriodOrbit(E,L)
    dt = halfPeriodUncertainty(E,L)
    return abs(dt/t)
end

##################################################
# Compute the orbit-averaged NR diffusion coefficients
# interpolate them beforehand
##################################################

nbrInt = 100

function averageDiffCoeffs(E::Float64, L::Float64, q::Float64, 
                           m_field::Float64)

    @assert (E>0.0 && L>0.0) "averageDiffCoeffs: E and L must be non-negative"

    halfperiod = halfPeriodOrbit(E,L)
    rmin, rmax = radiusBounds(E,L)

    rangerInt = range(rmin,length=nbrInt,rmax)
    tabrInt = collect(rangerInt)
    tabDiffCoeffsInt = zeros(Float64,5,nbrInt) # E,E2,L,L2,EL

    # Sample
    for indr=1:nbrInt
        rloc = tabrInt[indr]
        dEloc, dE2loc, dLloc, dL2loc, dEdLloc = localOrbitChange(rloc,E,L,q,m_field)
        tabDiffCoeffsInt[1][indr] = dEloc
        tabDiffCoeffsInt[2][indr] = dE2loc
        tabDiffCoeffsInt[3][indr] = dLloc
        tabDiffCoeffsInt[4][indr] = dL2loc
        tabDiffCoeffsInt[5][indr] = dEdLloc
    end

    # Interpolate
    intdEloc   = Interpolations.scale(interpolate(tabDiffCoeffsInt[1], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdE2loc  = Interpolations.scale(interpolate(tabDiffCoeffsInt[2], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdLloc   = Interpolations.scale(interpolate(tabDiffCoeffsInt[3], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdL2loc  = Interpolations.scale(interpolate(tabDiffCoeffsInt[4], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdEdLloc = Interpolations.scale(interpolate(tabDiffCoeffsInt[5], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)

    # Orbit-average
    dE   = quadgk(r->intdEloc(r)  /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)
    dE2  = quadgk(r->intdE2loc(r) /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)
    dL   = quadgk(r->intdLloc(r)  /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)
    dL2  = quadgk(r->intdL2loc(r) /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)
    dEdL = quadgk(r->intdEdLloc(r)/sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)

    dE   /= halfperiod
    dE2  /= halfperiod
    dL   /= halfperiod
    dL2  /= halfperiod
    dEdL /= halfperiod
    
    return dE, dE2, dL, dL2, dEdL
end