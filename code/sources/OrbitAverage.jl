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

# Returns "Unbounded orbit" if E<=0 and (rmin,rmax) if the parameters are allowed
# Check beforehands that the parameters are allowed, i.e!
# -> E,L are positive
# -> E < Ec(L)
function radiusBounds(E::Float64, L::Float64)
    if (E <= 0.0)
        return (0,"Unbounded orbit")
    else
        a = 4.0*E^2
        b = 4.0*(E^2-1.0+E*L^2)
        c = 4.0*E*L^2+L^4
        d = L^4
        nbRoot, rts = solveRealPoly3(a,b,c,d)

        @assert (nbRoot == 3) "radiusBounds: Forbidden parameters (E,L)"       
        rts = collect(rts)
        rts = sort(rts)
        rmin = sqrt(rts[2])
        rmax = sqrt(rts[3])
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

nbrInt = 10

function averageDiffCoeffs(E::Float64, L::Float64, q::Float64, 
                           m_field::Float64)

    @assert (E>0.0 && L>0.0) "averageDiffCoeffs: E and L must be non-negative"

    halfperiod = halfPeriodOrbit(E,L)
    rmin, rmax = radiusBounds(E,L)

    println(halfperiod)
    println(rmin)
    println(rmax)

    rangerInt = range(rmin,length=nbrInt,rmax)
    tabrInt = collect(rangerInt)
    tabDiffCoeffsInt = zeros(Float64,5,nbrInt) # E,E2,L,L2,EL

    println(tabrInt)

    # Sample
    for indr=1:nbrInt
        rloc = tabrInt[indr]
        dEloc, dE2loc, dLloc, dL2loc, dEdLloc = localOrbitChange(rloc,E,L,q,m_field)
        tabDiffCoeffsInt[1,indr] = dEloc
        tabDiffCoeffsInt[2,indr] = dE2loc
        tabDiffCoeffsInt[3,indr] = dLloc
        tabDiffCoeffsInt[4,indr] = dL2loc
        tabDiffCoeffsInt[5,indr] = dEdLloc
    end
 
    println(tabDiffCoeffsInt)
    
    # Interpolate
    intdEloc   = Interpolations.scale(interpolate(tabDiffCoeffsInt[1,:], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdE2loc  = Interpolations.scale(interpolate(tabDiffCoeffsInt[2,:], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdLloc   = Interpolations.scale(interpolate(tabDiffCoeffsInt[3,:], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdL2loc  = Interpolations.scale(interpolate(tabDiffCoeffsInt[4,:], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)
    intdEdLloc = Interpolations.scale(interpolate(tabDiffCoeffsInt[5,:], 
                                     BSpline(Cubic(Line(OnGrid())))),rangerInt)

    # Orbit-average
    dE   = quadgk(r->intdEloc(r)  /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
    dE2  = quadgk(r->intdE2loc(r) /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
    dL   = quadgk(r->intdLloc(r)  /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
    dL2  = quadgk(r->intdL2loc(r) /sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
    dEdL = quadgk(r->intdEdLloc(r)/sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]

    dE   /= halfperiod
    dE2  /= halfperiod
    dL   /= halfperiod
    dL2  /= halfperiod
    dEdL /= halfperiod
    
    return dE, dE2, dL, dL2, dEdL
end