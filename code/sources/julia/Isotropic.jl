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
    return DistFunction(E,0.1,0) # L=0.1 is hard-coded to get the correct if-condition
end

function localVelChangeIso!(r::Float64, v::Float64, m_field::Float64, 
                            PlummerTable::IntTable = PlummerTable_serial,
                            m_test::Float64=0.0)

    local E, cst, dvPar, dvPar2, dvTan2, int1, int3, K0, K1, K3
    let E, cst, dvPar, dvPar2, dvTan2, int1, int3, K0, K1, K3

    E = bindingEnergyIso(r,v)
    cst       = 16*pi^2*G^2*logCoulomb*m_field

    int1 = (Ea->velocityNorm(r,Ea)*DF_Iso(Ea))
    int3 = (Ea->velocityNorm(r,Ea)^3*DF_Iso(Ea))

    K0 = (2*E)^(9/2)/(21*pi^3)
    K1 = HCubature(int1,[E],[psi(r)],maxevals=MAXEVAL)[1]
    K3 = HCubature(int3,[E],[psi(r)],maxevals=MAXEVAL)[1]

    dvPar = -cst*(m_field+m_test)*K1
    dvPar2 = (2/3)*m_field*(K0+K3/v^3)
    dvTan2 = (2/3)*m_field*(2*K0+K1/v-K3/v^3)

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
