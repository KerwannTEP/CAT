##################################################
# Computation of the local velocity deflections
##################################################

using StaticArrays # To have access to static arrays
using Cuba

##################################################
# Actual Computation
##################################################

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

function _Ea(r::Float64 vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    return psi(r) - (1/2)*(v^2+vp^2-2*vp*(vr*cos(th)+vt*sin(th)*cos(phi)))
end

function _La(r::Float64 vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    return r*sqrt(vt^2+vp^2*sin(th)^2-2*vt*vp*sin(th)*cos(phi))
end

function dfadvr(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vt,vt,vp,th,phi)
    La = _La(r,vt,vt,vp,th,phi)
    return (-vr+vp*cos(th))*dDFdE(Ea,La,q)
end

function dfadvt(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vt,vt,vp,th,phi)
    La = _La(r,vt,vt,vp,th,phi)
    return (-vt+vp*sin(th)*cos(phi))*(dDFdE(Ea,La,q)-r/La *dDFdL(Ea,La,q))
end

function d2fadvr2(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                  th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vt,vt,vp,th,phi)
    La = _La(r,vt,vt,vp,th,phi)
    return -dDFdE(Ea,La,q) + (-vr+vp*cos(th))*d2DFdE2(Ea,La,q)
end

function d2fadvrdvt(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                    th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vt,vt,vp,th,phi)
    La = _La(r,vt,vt,vp,th,phi)
    return (-vr+vp*cos(th))*(-vt+vp*sin(th)*cos(phi))
           *(d2DFdE2(Ea,La,q) - r/La *d2DFdEdL(Ea,La,q))
end

function d2fadvt2(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                  th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vt,vt,vp,th,phi)
    La = _La(r,vt,vt,vp,th,phi)
    return -dDFdE(Ea,La,q) + r/La *-dDFdL (Ea,La,q) +(-vt+vp*sin(th)*cos(phi))^2
           *(d2DFdE2(Ea,La,q) - 2*r/La *d2DFdEdL(Ea,La,q)
             -r^2/La^3 *dDFdL(Ea,La,q) + r^2/La^2 *d2DFdL2(Ea,La,q))
end

function RosenbluthPotentials(r::Float64, vr::Float64, vt::Float64, q::Float64)
    #X = (v',x,phi) 
    vmax = _vmax(r,vr,vt)
    dhdvr     = 2*pi^2*vmax^2*cuhre((x, f) -> f[1] = x[1]  *sin(pi*x[2])    *dfadvr(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q),3,1)[1][1]
    dhdvt     = 2*pi^2*vmax^2*cuhre((x, f) -> f[1] = x[1]  *sin(pi*x[2])    *dfadvt(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q),3,1)[1][1]
    dgdvt     = 2*pi^2*vmax^4*cuhre((x, f) -> f[1] = x[1]^3*sin(pi*x[2])    *dfadvt(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q),3,1)[1][1]
    d2gdvr2   = 2*pi^2*vmax^4*cuhre((x, f) -> f[1] = x[1]^3*sin(pi*x[2])  *d2fadvr2(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q),3,1)[1][1]
    d2gdvrdvt = 2*pi^2*vmax^4*cuhre((x, f) -> f[1] = x[1]^3*sin(pi*x[2])*d2fadvrdvt(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q),3,1)[1][1]
    d2gdvt2   = 2*pi^2*vmax^4*cuhre((x, f) -> f[1] = x[1]^3*sin(pi*x[2])  *d2fadvt2(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q),3,1)[1][1]

    return dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2
end

function localVelChange(r::Float64, vr::Float64, vt::Float64,
                        q::Float64, m_field::Float64, m_test::Float64=0.0)

    dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2 = RosenbluthPotentials(r,vr,vt,q)
    cst = 4.0*pi*G^2*logCoulomb*m_field

    dvPar  = cst*(m_field+m_test) *((vr/v)*dhdvr+(vt/v)*dhdvt
    dvPar2 = cst* m_field         *((vr/v)^2*d2gdvr2+(2*vr*vt/v^2)*d2gdvrdvt
                                    +(vt/v)^2*d2gdvt2)
    dvTan2 = cst* m_field         *((vt/v)^2*d2gdvr2-(2*vr*vt/v^2)*d2gdvrdvt
                                    +(vr/v)^2*d2gdvt2+(1/vt)*dgdvt)

    return dvPar, dvPar2, dvTan2
end

function localOrbitChange(r::Float64, E::Float64, L::Float64, 
                          q::Float64, m_field::Float64, m_test::Float64=0.0)

    vr = radialVelocity(E,L,r)
    vt = tangentVelocity(E,L,r)
    v = sqrt(vr^2+vt^2)
    dvPar, dvPar2, dvTan2 = localVelChange(r,vr,vt,q,m_field,m_test)

    dE   = dvPar*(-v)  + dvPar2*(-1/2)    + dvTan2*(-1/2)
    dE2  =               dvPar2*(v^2)
    dL   = dvPar*(L/v)                    + dvTan2*(r^2/(4.0*L))
    dL2  =               dvPar2*(L^2/v^2) + dvTan2*(1/2 * (r^2 - L^2/v^2))
    dEdL =               dvPar2*(-L)
    return dE, dE2, dL, dL2, dEdL
end

function _ICuhre(r::Float64, vr::Float64, vt::Float64, q::Float64)
    #X = (v',x,phi) 
    vmax = _vmax(r,vr,vt)
    I1 = 2*pi^2*vmax*cuhre((x, f) -> f[1] = sin(pi*x[2])*dfPlummer(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q) ,3,1)[1][1]
    I2 = 2*pi^2*vmax^2*cuhre((x, f) -> f[1] = x[1]*sin(pi*x[2])*dfPlummer(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q) ,3,1)[1][1]

    return I1, I2
end

function localVelChangeOld(r::Float64, vr::Float64, vt::Float64,
                        q::Float64, m_field::Float64, m_test::Float64=0.0)

    cst = 4.0*pi*G^2*m_field
    dvPar  = cst*(m_field+m_test)*logCoulomb
    dvPar2 = cst* m_field 
    dvTan2 = cst* m_field        *(2.0*logCoulomb-1.0)
    I1, I2 = _ICuhre(r,vr,vt,q)
    dvPar  *= I1
    dvPar2 *= I2
    dvTan2 *= I2
    return dvPar, dvPar2, dvTan2
end

function localOrbitChangeOld(r::Float64, E::Float64, L::Float64, 
                          q::Float64, m_field::Float64, m_test::Float64=0.0)

    vr = radialVelocity(E,L,r)
    vt = tangentVelocity(E,L,r)
    v = sqrt(vr^2+vt^2)
    dvPar, dvPar2, dvTan2 = localVelChange(r,vr,vt,q,m_field,m_test)

    dE   = dvPar       + dvPar2*(-1/2)    + dvTan2*(-1/2)
    dE2  =               dvPar2*(v^2)
    dL   = dvPar*(-L/v)                   + dvTan2*(r^2/(4.0*L))
    dL2  =               dvPar2*(L^2/v^2) + dvTan2*(1/2 * (r^2 - L^2/v^2))
    dEdL =               dvPar2*(-L)
    return dE, dE2, dL, dL2, dEdL
end