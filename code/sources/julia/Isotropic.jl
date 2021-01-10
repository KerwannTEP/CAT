using StaticArrays # To have access to static arrays
using HCubature

MAXEVAL = 5000

function velocityNorm(r::Float64, E::Float64)
    return sqrt(2*(psi(r)-E))
end

function bindingEnergyIso(r::Float64, v::Float64)
    return psi(r)-(1/2)*v^2
end

function DF_Iso(E::Float64)
    return DistFunction(E,0.1,0.0) # L=0.1 is hard-coded to get the correct if-condition
end

function localVelChangeIso!(r::Float64, v::Float64, m_field::Float64, 
                            PlummerTable::IntTable = PlummerTable_serial,
                            m_test::Float64=0.0)

    local E, cst, dvPar, dvPar2, dvTan2, int1, int3, K0, K1, K3
    let E, cst, dvPar, dvPar2, dvTan2, int1, int3, K0, K1, K3

    E = bindingEnergyIso(r,v)
    cst = 16*pi^2*G^2*logCoulomb*m_field

    int1 = (Ea->velocityNorm(r,Ea[1])*DF_Iso(Ea[1]))
    int3 = (Ea->velocityNorm(r,Ea[1])^3*DF_Iso(Ea[1]))

    K0 = (2*E)^(9/2)/(21*pi^3)
    K1 = hcubature(int1,[E],[psi(r)],maxevals=MAXEVAL)[1]
    K3 = hcubature(int3,[E],[psi(r)],maxevals=MAXEVAL)[1]

    dvPar = -cst*(m_field+m_test)*K1
    dvPar2 = (2/3)*cst*m_field*(K0+K3/v^3)
    dvTan2 = (2/3)*cst*m_field*(2*K0+K1/v-K3/v^3)

    PlummerTable.dvPar[] = dvPar
    PlummerTable.dvPar2[] = dvPar2
    PlummerTable.dvTan2[] = dvTan2

    end
end

function localOrbitChangeIso!(r::Float64, E::Float64, L::Float64, m_field::Float64, 
                              PlummerTable::IntTable = PlummerTable_serial,
                              m_test::Float64=0.0)

    local v, dvPar, dvPar2, dvTan2, dE, dE2, dL, dL2, dEdL
    let v, dvPar, dvPar2, dvTan2, dE, dE2, dL, dL2, dEdL


    v = velocityNorm(r,E)

    localVelChangeIso!(r,v,m_field,PlummerTable,m_test)

    dvPar  = PlummerTable.dvPar[]
    dvPar2 = PlummerTable.dvPar2[]
    dvTan2 = PlummerTable.dvTan2[]

    dE   = dvPar*(-v)  + dvPar2*(-1/2)    + dvTan2*(-1/2)
    dE2  =               dvPar2*(v^2)
    dL   = dvPar*(L/v)                    + dvTan2*(r^2/(4.0*L))
    dL2  =               dvPar2*(L^2/v^2) + dvTan2*(1/2 * (r^2 - L^2/v^2))
    dEdL =               dvPar2*(-L)

    PlummerTable.dE[] = dE
    PlummerTable.dE2[] = dE2
    PlummerTable.dL[] = dL
    PlummerTable.dL2[] = dL2
    PlummerTable.dEdL[] = dEdL

    end
end


function averageDiffCoeffsIso!(E::Float64, L::Float64, m_field::Float64,
                               PlummerTable::IntTable = PlummerTable_serial)

    @assert (E>0.0 && L>0.0) "averageDiffCoeffs: E and L must be non-negative"

    IntTable_init!(PlummerTable)

    local halfperiod, rmin, rmax, rc, r1, r2, th1, th2, nbT1, nbT2, nbr,
          dE, dE2, dL, dL2, dEdL
    let halfperiod, rmin, rmax, rc, r1, r2, th1, th2, nbT1, nbT2, nbr,
        dE, dE2, dL, dL2, dEdL

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

        localOrbitChangeIso!(rp1,E,L,m_field,PlummerTable)

        dE   += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*PlummerTable.dE[] 
        dE2  += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*PlummerTable.dE2[] 
        dL   += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*PlummerTable.dL[]  
        dL2  += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*PlummerTable.dL2[]
        dEdL += dth1*cos(thp1)/abs(dpsiEffdr(rp1,L))*PlummerTable.dEdL[]
    end

    for iTh=1:nbT2
        thp2 = (iTh-0.5)*dth2
        rp2 = _orbitRadius(thp2,E,L,false) # max value (right of rc)

        halfperiod += dth2*cos(thp2)/abs(dpsiEffdr(rp2,L))
        
        localOrbitChangeIso!(rp2,E,L,m_field,PlummerTable)

        dE   += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*PlummerTable.dE[]  
        dE2  += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*PlummerTable.dE2[]  
        dL   += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*PlummerTable.dL[]  
        dL2  += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*PlummerTable.dL2[] 
        dEdL += dth1*cos(thp2)/abs(dpsiEffdr(rp2,L))*PlummerTable.dEdL[] 
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

        localOrbitChangeIso!(rp,E,L,m_field,PlummerTable)

        dE   += dr/sqrt(2*(psiEff(rp,L)-E))*PlummerTable.dE[] 
        dE2  += dr/sqrt(2*(psiEff(rp,L)-E))*PlummerTable.dE2[] 
        dL   += dr/sqrt(2*(psiEff(rp,L)-E))*PlummerTable.dL[] 
        dL2  += dr/sqrt(2*(psiEff(rp,L)-E))*PlummerTable.dL2[] 
        dEdL += dr/sqrt(2*(psiEff(rp,L)-E))*PlummerTable.dEdL[] 
    end
    dE   /= halfperiod
    dE2  /= halfperiod
    dL   /= halfperiod
    dL2  /= halfperiod
    dEdL /= halfperiod

    PlummerTable.dE[]   = dE
    PlummerTable.dE2[]  = dE2
    PlummerTable.dL[]   = dL
    PlummerTable.dL2[]  = dL2
    PlummerTable.dEdL[] = dEdL

    end
end
