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

nbT = 150
nbR = 150

function halfPeriodOrbit(E::Float64, L::Float64)
    rmin, rmax = radiusBounds(E,L)
    rc = maximizerPsiEff(L)
    r1, r2 = (rmin+rc)/2, (rc+rmax)/2
    th1, th2 = asin(sqrt(psiEff(r1,L)-E)), asin(sqrt(psiEff(r2,L)-E))
    dth1 = th1/nbT
    dth2 = th2/nbT
    dr = (r2-r1)/nbR
    halfperiod = 0.0
    for iTh=1:nbT
        thp1 = (iTh-0.5)*dth1
        thp2 = (iTh-0.5)*dth2
        rp1 = _orbitRadius(thp1,E,L,true) # min value (left of rc)
        rp2 = _orbitRadius(thp2,E,L,false) # max value (right of rc)
        halfPeriod += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))+dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))
    end
    halfPeriod *= sqrt(2)

    for iR=1:nbR
        rp = r1 + (iR-0.5)*dr
        halfPeriod += dr/sqrt(2*(psiEff(rp,L)-E))
    end
    return halfPeriod
end

##################################################
# Compute the orbit-averaged NR diffusion coefficients
# interpolate them beforehand
##################################################

# Transformation theta->r(theta)
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

nbrInt = 50

function averageDiffCoeffs(E::Float64, L::Float64, q::Float64, 
                           m_field::Float64)

    @assert (E>0.0 && L>0.0) "averageDiffCoeffs: E and L must be non-negative"

    halfperiod = 0.0
    rmin, rmax = radiusBounds(E,L)

    rangerInt = range(rmin,length=nbrInt,rmax)
    tabrInt = collect(rangerInt)
    tabDiffCoeffsInt = zeros(Float64,5,nbrInt) # E,E2,L,L2,EL
    
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

    rc = maximizerPsiEff(L)
    r1, r2 = (rmin+rc)/2, (rc+rmax)/2
    th1, th2 = asin(sqrt(psiEff(r1,L)-E)), asin(sqrt(psiEff(r2,L)-E))
    dth1 = th1/nbT
    dth2 = th2/nbT
    dr = (r2-r1)/nbR
    dE, dE2, dL, dL2, dEdL = 0.0, 0.0, 0.0, 0.0, 0.0

    # Orbit-average and computation of halfperiod
    for iTh=1:nbT
        thp1 = (iTh-0.5)*dth1
        thp2 = (iTh-0.5)*dth2
        rp1 = _orbitRadius(thp1,E,L,true)  # min value (left of rc)
        rp2 = _orbitRadius(thp2,E,L,false) # max value (right of rc)
        halfperiod += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))
        halfperiod += dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))
        dE   += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*intdEloc(rp1)   + dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))*intdEloc(rp2)
        dE2  += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*intdE2loc(rp1)  + dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))*intdE2loc(rp2)
        dL   += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*intdLloc(rp1)   + dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))*intdLloc(rp2)
        dL2  += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*intdL2loc(rp1)  + dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))*intdL2loc(rp2)
        dEdL += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*intdEdLloc(rp1) + dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))*intdEdLloc(rp2)
    end

    halfperiod *= sqrt(2)
    dE   *= sqrt(2)
    dE2  *= sqrt(2)
    dL   *= sqrt(2)
    dL2  *= sqrt(2)
    dEdL *= sqrt(2)

    for iR=1:nbR
        rp = r1 + (iR-0.5)*dr
        halfperiod += dr/sqrt(2*(psiEff(rp,L)-E))
        dE   += dr/sqrt(2*(psiEff(rp,L)-E))*intdEloc(rp)
        dE2  += dr/sqrt(2*(psiEff(rp,L)-E))*intdE2loc(rp)
        dL   += dr/sqrt(2*(psiEff(rp,L)-E))*intdLloc(rp)
        dL2  += dr/sqrt(2*(psiEff(rp,L)-E))*intdL2loc(rp)
        dEdL += dr/sqrt(2*(psiEff(rp,L)-E))*intdEdLloc(rp)
    end

    dE   /= halfperiod
    dE2  /= halfperiod
    dL   /= halfperiod
    dL2  /= halfperiod
    dEdL /= halfperiod
    
    return dE, dE2, dL, dL2, dEdL
end