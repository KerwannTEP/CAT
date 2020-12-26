##################################################
# Computation of the local velocity deflections
##################################################


# /!\ Impose a limit in Cubature


using HCubature
using StaticArrays # To have access to static arrays
using Cuba

nbPhi = 500
nbX = 500
nbResPoints = 100

##################################################
# Integration midpoints for the phi and x integrals
##################################################

tab_phi = SVector{nbPhi,Float64}([(pi/nbPhi)*(k-0.5) for k=1:nbPhi])
tab_x = SVector{nbX,Float64}([-1+(2/nbX)*(k-0.5) for k=1:nbX])

##################################################
# Actual Computation
##################################################

function _vmax(r::Float64, vr::Float64, vt::Float64)
    return vt + vr + sqrt(2.0*(vr*vt + psi(r)))
end

function dfSpherical(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                     x::Float64, phi::Float64, q::Float64)

    v  = sqrt(vr^2  + vt^2)
    Ea = psi(r) - 1/2*(v^2 + vp^2 - 2.0*vp*(vr*x + vt*sqrt(1.0-x^2)*cos(phi)))
    La = r*sqrt(vt^2 + vp^2*(1.0-x^2) - 2*vt*vp*sqrt(1.0-x^2)*cos(phi))
    DF = DistFunction(Ea, La, q)
    return DF
end

function _integralPhi(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                      x::Float64, q::Float64)
# integral over phi from 0 to 2pi
    intPhi = 0.0
    for k=1:nbX
 #       phik = tab_phi[k]
        fa = dfSpherical(r,vr,vt,vp,x,tab_phi[k],q) 
        intPhi += fa
    end
    intPhi *= 2*pi/nbPhi
    return intPhi
end

function _K(r::Float64, vr::Float64, vt::Float64, vp::Float64, q::Float64)
# integral over x from -1 to 1
    K = 0.0
    for k=1:nbX
 #       xk = tab_x[k]
        K += _integralPhi(r,vr,vt,vp,tab_x[k],q) 
    end
    K *= 2/nbX
    return K
end

function _KCubature(r::Float64, vr::Float64, vt::Float64, vp::Float64, q::Float64)
#integral over X = (x,phi) , -1<x<1 and 0<phi<2pi
    K = hcubature(x->dfSpherical(r,vr,vt,vp,x[1],x[2],q),[-1,0],[1,2*pi])
    return K[1]
end

function _I(r::Float64, vr::Float64, vt::Float64, q::Float64)
# integral I1 and _I2 of the notes
# I1 = int_0^vmax dv' K(v',r)
# I2 = int_0^vmax dv' v' K(v',r)
    vmax = _vmax(r,vr,vt)
    dvp = vmax/(nbResPoints)
    I1 = 0.0
    I2 = 0.0
    for k=1:nbResPoints
        K = _K(r,vr,vt,(k-0.5)*dvp,q)
        I1 += K
        I2 += K *(k-0.5)*dvp
    end
    return I1*dvp, I2*dvp
end

function DF_Cuba(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
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

function _ICubature(r::Float64, vr::Float64, vt::Float64, q::Float64)
    #X = (v',x,phi) 
    vmax = _vmax(r,vr,vt)
    I1 = hcubature(x->dfSpherical(r,vr,vt,x[1],x[2],x[3],q),[0,-1,0],[vmax,1,2*pi])
    I2 = hcubature(x->x[1]*dfSpherical(r,vr,vt,x[1],x[2],x[3],q),[0,-1,0],[vmax,1,2*pi])
    return I1[1], I2[1]
end

function _ICuba(r::Float64, vr::Float64, vt::Float64, q::Float64)
    #X = (v',x,phi) 
    vmax = _vmax(r,vr,vt)
    I1 = 2*pi^2*vmax*cuhre((x, f) -> f[1] = sin(pi*x[2])*DF_Cuba(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q) ,3,1)[1][1]
    I2 = 2*pi^2*vmax^2*cuhre((x, f) -> f[1] = x[1]*sin(pi*x[2])*DF_Cuba(r,vr,vt,vmax*x[1],pi*x[2],2*pi*x[3],q) ,3,1)[1][1]

    return I1, I2
end

function localVelChange(r::Float64, vr::Float64, vt::Float64,
                        q::Float64, m_field::Float64, m_test::Float64=0.0)

    cst = 4.0*pi*G^2*m_field
    dvPar  = cst*(m_field+m_test)*logCoulomb
    dvPar2 = cst* m_field 
    dvTan2 = cst* m_field        *(2.0*logCoulomb-1.0)
    I1, I2 = _ICuba(r,vr,vt,q)
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