##################################################
# Computation of the local velocity deflections
##################################################

"""
    _vmax(r,vr,vt,sinth,costh,cosph)

Computes the maximum allowed velocity in Rosenbluth integrals, for a test star at position `r` and velocity `(vr,vt)`.
At fixed integration variable (over field) `(θ,φ)` such that `sin(θ)=sinth`, `cos(θ)=costh` and `cos(φ)=cosph`.

# Remarks
- Integrands in Rosenbluth integrals vanishes for v > _vmax.

# Arguments
- `r::Float64`: Radial position of the test star.
- `vr::Float64`: Radial velocity of the test star.
- `vt::Float64`: Tangential velocity of the test star.
- `sinth::Float64`: Value `sin(θ)=sinth` corresponding to `θ` integration variable.
- `costh::Float64`: Value `cos(θ)=costh` corresponding to `θ` integration variable.
- `cosph::Float64`: Value `cos(φ)=cosph` corresponding to `φ` integration variable.
"""
function _vmax(r::Float64, vr::Float64, vt::Float64,
    sinth::Float64,costh::Float64,cosph::Float64)

    vSq = vr^2 + vt^2
    deltav = (vr*costh+vt*sinth*cosph)^2 - (2*psi(r)+vSq)


    if (deltav > 0.0)
        return (vr*costh+vt*sinth*cosph) + sqrt(deltav)
    else
        return 0.0
    end

end

"""
    _EaLa(r,var2,vat2)

Computes the energy and angular momentum used in Rosenbluth integrands.


# Arguments
- `r::Float64`: Radial position of the test star.
- `var2::Float64`: Square radial velocity in Rosenbluth integrand.
- `vat2::Float64`: Square tangential velocity in Rosenbluth integrand.
"""
function _EaLa(r::Float64, var2::Float64, vat2::Float64)
    Ea = psi(r) + (1/2)*(var2+vat2)
    La = r*sqrt(vat2)
    return Ea, La
end

"""
    _fa(r,vr,vt,vp,snth,csth,snph,csph)

Distribution function used in Rosenbluth integrand


# Arguments
- `r::Float64`: Radial position of the test star.
- `vr::Float64`: Radial velocity of the test star.
- `vt::Float64`: Tangential velocity of the test star.
- `vp::Float64`: Velocity integration variable.
- `snth::Float64`: Value `sin(θ)=snth` corresponding to `θ` integration variable.
- `csth::Float64`: Value `cos(θ)=csth` corresponding to `θ` integration variable.
- `snph::Float64`: Value `sin(φ)=snph` corresponding to `φ` integration variable.
- `csph::Float64`: Value `cos(φ)=csph` corresponding to `φ` integration variable.
"""
function _fa(r::Float64, vr::Float64, vt::Float64, vp::Float64,
                snth::Float64,csth::Float64,snph::Float64,csph::Float64)

    pref_r = vr-vp*csth
    pref_t = vt-vp*snth*csph
    pref_t2 = vp*snth*snph

    var2 = pref_r^2
    vat2 = pref_t^2+pref_t2^2
    Ea, La = _EaLa(r,var2,vat2)
    DF = _F(Ea,La)

    return DF
end


"""
    RosenbluthPotentials(r,vr,vt,[nbK=nbK_default])

Computation of the Rosenbluth potential and its derivatives.
Returns a tuple `(dh/dvr, dh/dvt, dg/dvt, d2g/dvr2, d2g/dvrdvt, d2g/dvt2)`.

# Arguments
- `r::Float64`: Radial position of the test star.
- `vr::Float64`: Radial velocity of the test star.
- `vt::Float64`: Tangential velocity of the test star.
- `nbK::Int64`: Sampling number of Rosenbluth integrals.
"""
function RosenbluthPotentials(r::Float64, vr::Float64, vt::Float64, nbK::Int64=nbK_default)
    #X = (v',x,phi)

    sumhr = 0
    sumht = 0
    sumgt = 0
    sumh = 0
    sumgrrInt = 0
    sumgrt = 0
    sumgttInt = 0
    for ith=1:nbK
        thtemp = (2*ith-1)/(2*nbK)
        sinthtemp,costhtemp = sincos(PI*thtemp)
        sumhr_ph = 0
        sumht_ph = 0
        sumgt_ph = 0
        sumh_ph = 0
        sumgrrInt_ph = 0
        sumgrt_ph = 0
        sumgttInt_ph = 0

        for iph=1:nbK
            phtemp = (2*iph-1)/(2*nbK)
            sinphtemp, cosphtemp = sincos(2*PI*phtemp)
            vmax = _vmax(r,vr,vt,sinthtemp,costhtemp,cosphtemp)

            sumhr_v = 0
            sumht_v = 0
            sumgt_v = 0
            sumh_v = 0
            sumgrrInt_v = 0
            sumgrt_v = 0
            sumgttInt_v = 0

            for iv=1:nbK
                vtemp = (2*iv-1)/(2*nbK)
                floc = _fa(r,vr,vt,vmax*vtemp,sinthtemp,costhtemp,sinphtemp,cosphtemp)
                sumhr_v += -floc
                sumht_v += -floc
                sumgt_v += vtemp^2*floc

                # compute h and the other integrals separately
                sumh_v += vtemp*floc
                sumgrrInt_v += -vtemp*floc
                sumgrt_v += -vtemp*floc
                sumgttInt_v += -vtemp*floc
            end
            sumhr_v *= vmax
            sumht_v *= vmax
            sumgt_v *= vmax^3
            sumh_v *= vmax^2
            sumgrrInt_v *= vmax^2
            sumgrt_v *= vmax^2
            sumgttInt_v *= vmax^2

            sumhr_ph += sumhr_v
            sumht_ph += sumht_v*cosphtemp
            sumgt_ph += sumgt_v*cosphtemp
            sumh_ph += sumh_v
            sumgrrInt_ph += sumgrrInt_v
            sumgrt_ph += sumgrt_v*cosphtemp
            sumgttInt_ph += sumgttInt_v*cosphtemp^2
        end
        sumhr += sinthtemp*costhtemp*sumhr_ph
        sumht += sinthtemp^2*sumht_ph
        sumgt += sinthtemp^2*sumgt_ph
        sumh += sinthtemp*sumh_ph
        sumgrrInt += sinthtemp*costhtemp^2*sumgrrInt_ph
        sumgrt += sinthtemp^2*costhtemp*sumgrt_ph
        sumgttInt += sinthtemp^3*sumgttInt_ph
    end
    sumhr *= PI*2*PI/(nbK^3)
    sumht *= PI*2*PI/(nbK^3)
    sumgt *= PI*2*PI/(nbK^3)
    sumh *= PI*2*PI/(nbK^3)
    sumgrrInt *= PI*2*PI/(nbK^3)
    sumgrt *= PI*2*PI/(nbK^3)
    sumgttInt *= PI*2*PI/(nbK^3)

    return sumhr, sumht, sumgt, sumh+sumgrrInt, sumgrt, sumh+sumgttInt
end

"""
    localVelChange(r,vr,vt,m_field,[nbK=nbK_default,m_test=m_field])

Computation of the local velocity deflections along and perpendicular to the trajectory of the test star.
Returns a tuple `(dvPar, dvPar2, dvPerp2)`.

# Arguments
- `r::Float64`: Radial position of the test star.
- `vr::Float64`: Radial velocity of the test star.
- `vt::Float64`: Tangential velocity of the test star.
- `m_field::Float64`: Mass of field star.
- `nbK::Int64`: Sampling number of Rosenbluth integrals.
- `m_test::Float64`: Mass of test star.
"""
function localVelChange(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64, nbK::Int64=nbK_default,
                        m_test::Float64=m_field)

    dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2 = RosenbluthPotentials(r,vr,vt,nbK)

    cst       = 4.0*PI*_G^2*logCoulomb
    v         = sqrt(vr^2+vt^2)

    if (vt>0)
        vr_v   = vr/v
        vt_v   = vt/v

        dvPar   = cst*(m_field+m_test) *(vr_v*dhdvr+vt_v*dhdvt)
        dvPar2  = cst* m_field         *(vr_v^2*d2gdvr2+(2*vr_v*vt_v)*d2gdvrdvt
                                        +vt_v^2*d2gdvt2)
        dvPerp2 = cst* m_field         *(vt_v^2*d2gdvr2-(2*vr_v*vt_v)*d2gdvrdvt
                                        +vr_v^2*d2gdvt2+(1/vt)*dgdvt)
    else # vt = 0

        dvPar   = cst*(m_field+m_test) *dhdvr
        dvPar2  = cst* m_field         *d2gdvr2
        dvPerp2 = cst* m_field         *2*d2gdvt2

    end

    return dvPar, dvPar2, dvPerp2
end

"""
    localOrbitChange(r,E,L,m_field,[nbK=nbK_default,m_test=m_field])

Computation of the local energy and angular momentum deflections.
Returns a tuple `(dE, dE2, dL, dEdL, dL2)`.

# Arguments
- `r::Float64`: Radial position of the test star.
- `E::Float64`: Energy of the test star.
- `L::Float64`: Angular momentum of the test star.
- `m_field::Float64`: Mass of field star.
- `nbK::Int64`: Sampling number of Rosenbluth integrals.
- `m_test::Float64`: Mass of test star.
"""
function localOrbitChange(r::Float64, E::Float64, L::Float64,
                           m_field::Float64, nbK::Int64=nbK_default,
                           m_test::Float64=m_field)

    vr = radialVelocity(r,E,L)
    vt = tangentVelocity(r,E,L)
    v = sqrt(vr^2+vt^2)

    dvPar, dvPar2, dvPerp2 = localVelChange(r,vr,vt,m_field,nbK,m_test)

    dE   = dvPar*(v)   + dvPar2*(1/2)     + dvPerp2*(1/2)
    dE2  =               dvPar2*(v^2)
    dL   = dvPar*(L/v)                    + dvPerp2*(r^2/(4.0*L))
    dL2  =               dvPar2*(L^2/v^2) + dvPerp2*(1/2 * (r^2 - L^2/v^2))
    dEdL =               dvPar2*(L)

    return dE, dE2, dL, dEdL, dL2
end
