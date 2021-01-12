##################################################
# Computation of the local velocity deflections
##################################################

using StaticArrays # To have access to static arrays
using HCubature

##################################################
# Actual Computation
##################################################

PlummerTable_serial = IntTable_create!()

function _vmax(r::Float64, vr::Float64, vt::Float64)
    return vt + vr + sqrt(2.0*(vr*vt + psi(r)))
end

function dfPlummer(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                 th::Float64, phi::Float64, q::Float64)
    v  = sqrt(vr^2  + vt^2)
    csth = cos(th)
    snth = sin(th)
    csph = cos(phi)
    Ea = psi(r) - 1/2*(v^2 + vp^2 - 2.0*vp*(vr*csth + vt*snth*csph))
    La = r*sqrt(vt^2 + vp^2*snth^2 - 2*vt*vp*snth*csph)
    DF = DistFunction(Ea, La, q)
    return DF
end

function _Ea(r::Float64, vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    v  = sqrt(vr^2  + vt^2)
    return psi(r) - (1/2)*(v^2+vp^2-2*vp*(vr*cos(th)+vt*sin(th)*cos(phi)))
end

function _La(r::Float64, vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    return r*sqrt(vt^2+vp^2*sin(th)^2-2*vt*vp*sin(th)*cos(phi))
end

function dfadvr(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vr+vp*cos(th))*dDFdE(Ea,La,q)
end

function dfadvt(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vt+vp*sin(th)*cos(phi))*(dDFdE(Ea,La,q)-r/La *dDFdL(Ea,La,q))
end

function d2fadvr2(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                  th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return -dDFdE(Ea,La,q) + (-vr+vp*cos(th))*d2DFdE2(Ea,La,q)
end

function d2fadvrdvt(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                    th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vr+vp*cos(th))*(-vt+vp*sin(th)*cos(phi))
           *(d2DFdE2(Ea,La,q) - r/La *d2DFdEdL(Ea,La,q))
end

function d2fadvt2(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                  th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return -dDFdE(Ea,La,q) + r/La *-dDFdL(Ea,La,q) +(-vt+vp*sin(th)*cos(phi))^2
           *(d2DFdE2(Ea,La,q) - 2*r/La *d2DFdEdL(Ea,La,q)
             -r^2/La^3 *dDFdL(Ea,La,q) + r^2/La^2 *d2DFdL2(Ea,La,q))
end

MAXEVAL = 5000


function RosenbluthPotentials!(r::Float64, vr::Float64, vt::Float64, q::Float64,
                               PlummerTable::IntTable = PlummerTable_serial)
    #X = (v',x,phi) 
    vmax = _vmax(r,vr,vt)

    local int_hr, int_ht, int_gt, int_grr, int_grt, int_gtt
    let int_hr, int_ht, int_gt, int_grr, int_grt, int_gtt

    int_hr = (x->x[1]  *sin(pi*x[2]) *dfadvr(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q))
    int_ht = (x-> x[1]  *sin(pi*x[2])    *dfadvt(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q))
    int_gt = (x-> x[1]^3*sin(pi*x[2])    *dfadvt(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q))
    int_grr = (x->x[1]^3*sin(pi*x[2])  *d2fadvr2(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q))
    int_grt = (x->x[1]^3*sin(pi*x[2])*d2fadvrdvt(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q))
    int_gtt = (x->x[1]^3*sin(pi*x[2])  *d2fadvt2(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q))


    PlummerTable.dhdvr[]     = 2*pi^2*vmax^2*hcubature(int_hr,[0,0,0],[1,1,1],maxevals=MAXEVAL)[1]
    PlummerTable.dhdvt[]     = 2*pi^2*vmax^2*hcubature(int_ht,[0,0,0],[1,1,1],maxevals=MAXEVAL)[1]
    PlummerTable.dgdvt[]     = 2*pi^2*vmax^4*hcubature(int_gt,[0,0,0],[1,1,1],maxevals=MAXEVAL)[1]
    PlummerTable.d2gdvr2[]   = 2*pi^2*vmax^4*hcubature(int_grr,[0,0,0],[1,1,1],maxevals=MAXEVAL)[1]
    PlummerTable.d2gdvrdvt[] = 2*pi^2*vmax^4*hcubature(int_grt,[0,0,0],[1,1,1],maxevals=MAXEVAL)[1]
    PlummerTable.d2gdvt2[]   = 2*pi^2*vmax^4*hcubature(int_gtt,[0,0,0],[1,1,1],maxevals=MAXEVAL)[1]

    end

end

function localVelChange!(r::Float64, vr::Float64, vt::Float64,
                        q::Float64, m_field::Float64, PlummerTable::IntTable = PlummerTable_serial,
                        m_test::Float64=0.0)


    RosenbluthPotentials!(r,vr,vt,q,PlummerTable)    

    local v, cst, dvPar, dvPar2, dvTan2, dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2
    let v, cst, dvPar, dvPar2, dvTan2, dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2

    cst       = 4.0*pi*G^2*logCoulomb*m_field
    v         = sqrt(vr^2+vt^2)
    dhdvr     = PlummerTable.dhdvr[]
    dhdvt     = PlummerTable.dhdvt[]
    dgdvt     = PlummerTable.dgdvt[] 
    d2gdvr2   = PlummerTable.d2gdvr2[]
    d2gdvrdvt = PlummerTable.d2gdvrdvt[]
    d2gdvt2   = PlummerTable.d2gdvt2[]

    dvPar  = -cst*(m_field+m_test) *((vr/v)*dhdvr+(vt/v)*dhdvt)
    dvPar2 = cst* m_field         *((vr/v)^2*d2gdvr2+(2*vr*vt/v^2)*d2gdvrdvt
                                    +(vt/v)^2*d2gdvt2)
    dvTan2 = cst* m_field         *((vt/v)^2*d2gdvr2-(2*vr*vt/v^2)*d2gdvrdvt
                                    +(vr/v)^2*d2gdvt2+(1/vt)*dgdvt)

    PlummerTable.dvPar[] = dvPar
    PlummerTable.dvPar2[] = dvPar2
    PlummerTable.dvTan2[] = dvTan2

    end

end

function localOrbitChange!(r::Float64, E::Float64, L::Float64,
                           q::Float64, m_field::Float64, PlummerTable::IntTable = PlummerTable_serial,
                           m_test::Float64=0.0)

    local v, vr, vt, dvPar, dvPar2, dvTan2, dE, dE2, dL, dL2, dEdL
    let v, vr, vt, dvPar, dvPar2, dvTan2, dE, dE2, dL, dL2, dEdL

    vr = radialVelocity(r,E,L)
    vt = tangentVelocity(r,E,L)
    v = sqrt(vr^2+vt^2)

    
    localVelChange!(r,vr,vt,q,m_field,PlummerTable,m_test)

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
