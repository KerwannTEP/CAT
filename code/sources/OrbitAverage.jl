##################################################
# Orbit-averagering
##################################################

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

nbT = 100
nbR = 100

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
        halfperiod += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))+dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))
    end
    halfperiod *= sqrt(2)

    for iR=1:nbR
        rp = r1 + (iR-0.5)*dr
        halfperiod += dr/sqrt(2*(psiEff(rp,L)-E))
    end
    return halfperiod
end

##################################################
# Compute the orbit-averaged NR diffusion coefficients
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

function averageDiffCoeffs(E::Float64, L::Float64, q::Float64, 
                           m_field::Float64)

    @assert (E>0.0 && L>0.0) "averageDiffCoeffs: E and L must be non-negative"

    halfperiod = 0.0
    rmin, rmax = radiusBounds(E,L)
    rc = maximizerPsiEff(L)
    r1, r2 = (rmin+rc)/2, (rc+rmax)/2
    th1, th2 = asin(sqrt(psiEff(r1,L)-E)), asin(sqrt(psiEff(r2,L)-E))   

    nbT1 = 20
    nbT2 = 20
    nbr = 30

    dth1 = th1/nbT1
    dth2 = th2/nbT2
    dr = (r2-r1)/nbr
    dE, dE2, dL, dL2, dEdL = 0.0, 0.0, 0.0, 0.0, 0.0

    # Orbit-average and computation of halfperiod
    for iTh=1:nbT1
        thp1 = (iTh-0.5)*dth1
        rp1 = _orbitRadius(thp1,E,L,true)  # min value (left of rc)
        halfperiod += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))

        dEloc, dE2loc, dLloc, dL2loc, dEdLloc = localOrbitChange(rp1,E,L,q,m_field)

        dE   += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*dEloc 
        dE2  += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*dE2loc 
        dL   += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*dLloc  
        dL2  += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*dL2loc
        dEdL += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*dEdLloc
    end

    for iTh=1:nbT2
        thp2 = (iTh-0.5)*dth2
        rp2 = _orbitRadius(thp2,E,L,false) # max value (right of rc)
        halfperiod += dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))
        dEloc, dE2loc, dLloc, dL2loc, dEdLloc = localOrbitChange(rp2,E,L,q,m_field)

        dE   += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*dEloc 
        dE2  += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*dE2loc 
        dL   += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*dLloc  
        dL2  += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*dL2loc
        dEdL += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*dEdLloc
    end
    halfperiod *= sqrt(2)
    dE   *= sqrt(2)
    dE2  *= sqrt(2)
    dL   *= sqrt(2)
    dL2  *= sqrt(2)
    dEdL *= sqrt(2)

    for iR=1:nbr
        rp = r1 + (iR-0.5)*dr
        halfperiod += dr/sqrt(2*(psiEff(rp,L)-E))
        dEloc, dE2loc, dLloc, dL2loc, dEdLloc = localOrbitChange(rp,E,L,q,m_field)

        dE   += dr/sqrt(2*(psiEff(rp,L)-E))*dEloc
        dE2  += dr/sqrt(2*(psiEff(rp,L)-E))*dE2loc
        dL   += dr/sqrt(2*(psiEff(rp,L)-E))*dLloc
        dL2  += dr/sqrt(2*(psiEff(rp,L)-E))*dL2loc
        dEdL += dr/sqrt(2*(psiEff(rp,L)-E))*dEdLloc
    end
    dE   /= halfperiod
    dE2  /= halfperiod
    dL   /= halfperiod
    dL2  /= halfperiod
    dEdL /= halfperiod
    return dE, dE2, dL, dL2, dEdL
end