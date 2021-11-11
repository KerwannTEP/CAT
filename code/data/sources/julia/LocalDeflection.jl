##################################################
# Computation of the local velocity deflections
##################################################

function _deltav(r::Float64, vr::Float64, vt::Float64,
    th::Float64, phi::Float64)

    vSq = vr^2 + vt^2
    return (vr*cos(th)+vt*sin(th)*cos(phi))^2 - (2*psi(r)+vSq)
end

function _vmax(r::Float64, vr::Float64, vt::Float64,
    th::Float64, phi::Float64)

    deltav = _deltav(r,vr,vt,th,phi)

    if (deltav > 0.0)
        return (vr*cos(th)+vt*sin(th)*cos(phi)) + sqrt(deltav)
    else
        return 0.0
    end
end

function _vmaxNew(r::Float64, vr::Float64, vt::Float64,
    sinth::Float64,costh::Float64,cosph::Float64)

    vSq = vr^2 + vt^2
    deltav = (vr*costh+vt*sinth*cosph)^2 - (2*psi(r)+vSq)


    if (deltav > 0.0)
        return (vr*costh+vt*sinth*cosph) + sqrt(deltav)
    else
        return 0.0
    end

end

function _EaLa(r::Float64, var2::Float64, vat2::Float64)
    Ea = psi(r) + (1/2)*(var2+vat2)
    La = r*sqrt(vat2)
    return Ea, La
end



function _fa(r::Float64, vr::Float64, vt::Float64, vp::Float64,
                snth::Float64,csth::Float64,snph::Float64,csph::Float64)

    # snth, csth = sincos(th)
    # snph, csph = sincos(phi)
    pref_r = vr-vp*csth
    pref_t = vt-vp*snth*csph
    pref_t2 = vp*snth*snph

    var2 = pref_r^2
    vat2 = pref_t^2+pref_t2^2

    Ea, La = _EaLa(r,var2,vat2)

    DF = _F(Ea,La)

    return DF
end



function RosenbluthPotentialsSum(r::Float64, vr::Float64, vt::Float64, nbK::Int=nbK_default)
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

            vmax = _vmaxNew(r,vr,vt,sinthtemp,costhtemp,cosphtemp)

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

function localVelChange(r::Float64, vr::Float64, vt::Float64,
                        m_field::Float64, nbK::Int=nbK_default,
                        m_test::Float64=m_field)

    dhdvr, dhdvt, dgdvt, d2gdvr2, d2gdvrdvt, d2gdvt2 = RosenbluthPotentialsSum(r,vr,vt,nbK)

    cst       = 4.0*PI*_G^2*logCoulomb
    v         = sqrt(vr^2+vt^2)

    if (vt>0)
        vr_v   = vr/v
        vt_v   = vt/v

        dvPar  = cst*(m_field+m_test) *(vr_v*dhdvr+vt_v*dhdvt)
        dvPar2 = cst* m_field         *(vr_v^2*d2gdvr2+(2*vr_v*vt_v)*d2gdvrdvt
                                        +vt_v^2*d2gdvt2)
        dvPerp2 = cst* m_field         *(vt_v^2*d2gdvr2-(2*vr_v*vt_v)*d2gdvrdvt
                                        +vr_v^2*d2gdvt2+(1/vt)*dgdvt)
    else # vt = 0

        dvPar  = cst*(m_field+m_test) *dhdvr
        dvPar2 = cst* m_field         *d2gdvr2
        dvPerp2 = cst* m_field         *2*d2gdvt2

    end

    return dvPar, dvPar2, dvPerp2
end


function localOrbitChange(r::Float64, E::Float64, L::Float64,
                           m_field::Float64, nbK::Int=nbK_default,
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
