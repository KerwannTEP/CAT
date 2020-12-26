##################################################
# Orbit-averagering
##################################################

using QuadGK
using Interpolations

##################################################
# Determine the nature of the orbit
# Use positive binding energy
##################################################

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
        solver = solveRealPoly3(a,b,c,d)

        @assert (solver[1] == 3) "radiusBounds: Forbidden parameters (E,L)" 
        rt1, rt2, rt3 = solver[2], solver[3], solver[4]
        if (rt1 < 0.0)
            rmin = rt2
            rmax = rt3
        else
            rmin = rt1
            if (rt2 < 0.0)
                rmax = rt3
            else
                rmax = rt2
            end
        end
        if (rmax < rmin)
            rmin, rmax = rmax, rmin
        end
        rmin = sqrt(rmin)
        rmax = sqrt(rmax)
        return rmin, rmax
    end
end

##################################################
# Compute the period of an orbit
##################################################

function _orbitRadius(th::Float64, E::Float64, L::Float64, left::Bool)
    a = 4.0*(E+sin(th)^2)^2
    b = 4.0*((E+sin(th)^2)^2-1.0+(E+sin(th)^2)*L^2)
    c = 4.0*(E+sin(th)^2)*L^2+L^4
    d = L^4
    solver = solveRealPoly3(a,b,c,d)
    @assert (solver[1] == 3) "_orbitRadius: Forbidden angle theta" 
    rt1, rt2, rt3 = solver[2], solver[3], solver[4]
    if (rt1 < 0.0)
        rmin = rt2
        rmax = rt3
    else
        rmin = rt1
        if (rt2 < 0.0)
            rmax = rt3
        else
            rmax = rt2
        end
    end
    if (rmax < rmin)
        rmin, rmax = rmax, rmin
    end
    rmin = sqrt(rmin)
    rmax = sqrt(rmax)
    if (left)
        return rmin
    else
        return rmax
    end
end


nbT = 100
nbR = 100

function testChangeVar(E::Float64, L::Float64)
    rmin, rmax = radiusBounds(E,L)
    rc = maximizerPsiEff(L)
    r1, r2 = (rmin+rc)/2, (rc+rmax)/2
    th1, th2 = asin(sqrt(psiEff(r1,L)-E)), asin(sqrt(psiEff(r2,L)-E))
    dth1 = th1/nbT
    dth2 = th2/nbT
    dr = (r2-r1)/nbR
    t1, t2, t3 = 0.0, 0.0, 0.0
    for iTh=1:nbT
        thp1 = (iTh-0.5)*dth1
        thp2 = (iTh-0.5)*dth2
       rp1 = _orbitRadius(thp1,E,L,true) # min value (left of rc)
       rp2 = _orbitRadius(thp2,E,L,false) # max value (right of rc)
        t1 += cos(thp1)/abs(dpsiEffdr(rp1,L))
        t2 += cos(thp2)/abs(dpsiEffdr(rp2,L))
    end
    t1 *= dth1*sqrt(2)
    t2 *= dth2*sqrt(2)

    for iR=1:nbR
        rp = r1 + (iR-0.5)*dr
        t3 += 1/sqrt(2*(psiEff(rp,L)-E))
    end
    t3 *= dr
  #  t3 = quadgk(r->1/sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
 #   println(t1)
  #  println(t2)
   # println(t3)
    return t1+t2+t3
end

function halfPeriodOrbit(E::Float64, L::Float64)
    rmin, rmax = radiusBounds(E,L)
    t3 = quadgk(r->1/sqrt(2.0*abs(psiEff(r,L)-E)),rmin,rmax)[1]
    return t3
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

nbrInt = 50

function averageDiffCoeffs(E::Float64, L::Float64, q::Float64, 
                           m_field::Float64)

    @assert (E>0.0 && L>0.0) "averageDiffCoeffs: E and L must be non-negative"

    halfperiod = halfPeriodOrbit(E,L)
    rmin, rmax = radiusBounds(E,L)

    rangerInt = range(rmin,length=nbrInt,rmax)
    tabrInt = collect(rangerInt)
    tabDiffCoeffsInt = zeros(Float64,5,nbrInt) # E,E2,L,L2,EL
    
#    dEloc, dE2loc, dLloc, dL2loc, dEdLloc = 0.0, 0.0, 0.0, 0.0, 0.0

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