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


function _ICuhre(r::Float64, vr::Float64, vt::Float64, q::Float64)
    #X = (v',x,phi) 
    vmax = _vmax(r,vr,vt)
    I1 = 2*pi^2*vmax*cuhre((x, f) -> f[1] = sin(pi*x[2])*dfPlummer(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q) ,3,1)[1][1]
    I2 = 2*pi^2*vmax^2*cuhre((x, f) -> f[1] = x[1]*sin(pi*x[2])*dfPlummer(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q) ,3,1)[1][1]

    return I1, I2
end

function localVelChange(r::Float64, vr::Float64, vt::Float64,
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

function localOrbitChange(r::Float64, E::Float64, L::Float64, 
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