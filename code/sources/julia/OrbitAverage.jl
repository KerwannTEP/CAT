##################################################
# Orbit-averagering of the diffusion coefficients
# In (E,L) space
##################################################
using PolynomialRoots


function radius_s_bounds(E::Float64, L::Float64)
    if (E >= 0.0)
        return "Unbounded orbit"
    elseif (L > Lc(_tE(E)))
        return "Not a possible orbit"
    else
        tE = _tE(E)
        tL = _tL(L)
        if (L != 0.0)
            rts = sort(real(roots([1,-(tE-tL^2/2),-1,tE])))
            return Float64(rts[2]), Float64(rts[3])
        else
            return 1.0, 1/tE
        end
    end
end


# does not check if orbit is correct
# only use when sure the orbit is correct
function radius_s_bounds_unsafe(E::Float64, L::Float64)

    tE = _tE(E)
    tL = _tL(L)
    if (L != 0.0)
        rts = sort(real(roots([1,-(tE-tL^2/2),-1,tE])))
        return Float64(rts[2]), Float64(rts[3])
    else
        return 1.0, 1/tE
    end

end


function sma_eff(sp::Float64, sa::Float64)
    return (sa+sp)/2
end

function ecc_eff(sp::Float64, sa::Float64)
    return (sa-sp)/(sa+sp)
end

#tr=r/b
function tr_eff_anomaly(u::Float64, a::Float64, ecc::Float64)
    s = a*(1+ecc*u*(1.5-u^2/2))

    return sqrt(abs(s^2-1)) # use absolute value for when s->1


end

function r_eff_anomaly(u::Float64, a::Float64, ecc::Float64)
    return _b*tr_eff_anomaly(u,a,ecc)
end

function dth1du(u::Float64, sp::Float64, sa::Float64)
    num = sa*(-2+u)*(1+u)^2-sp*(2-3*u+u^3)
    denumSq = ((sa^2*sp*(-2+u)*(1+u)^2-sp*(6-3*u+u^3)
               +sa*(-6-3*u+u^3-sp^2*(2-3*u+u^3)))/
               (sa*sp*(sa+sp)*(sa*(-2+u)*(1+u)^2-sp*(2-3*u+u^3))))
    return -3/(4*sqrt(2))*1/sqrt(4-u^2)*num/sqrt(denumSq)
end

function halfperiod(E::Float64, L::Float64, nbAvr::Int64=nbAvr_default)
    sp, sa = radius_s_bounds(E,L)
    halfperiod = 0.0

    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        jac_loc = dth1du(uloc,sp,sa)
        halfperiod += jac_loc
    end

    halfperiod *= (1/_Omega0)*(2/nbAvr)

    return halfperiod
end

function avr_orbit_coefficients(E::Float64, L::Float64, m_field::Float64,
                                nbK::Int=nbK_default, nbAvr::Int64=nbAvr_default,
                                m_test::Float64=m_field)
    sp, sa = radius_s_bounds(E,L)
    a = sma_eff(sp,sa)
    ecc = ecc_eff(sp,sa)

    avrDE = 0.0
    avrDL = 0.0
    avrDEE = 0.0
    avrDEL = 0.0
    avrDLL = 0.0
    halfperiod = 0.0

    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        rloc = r_eff_anomaly(uloc,a,ecc)
        jac_loc = dth1du(uloc,sp,sa)
        dEloc, dE2loc, dLloc, dEdLloc, dL2loc = localOrbitChangeReg(rloc,E,L,m_field,nbK,m_test)

        avrDE += jac_loc*dEloc
        avrDL += jac_loc*dLloc
        avrDEE += jac_loc*dE2loc
        avrDEL += jac_loc*dEdLloc
        avrDLL += jac_loc*dL2loc
        halfperiod += jac_loc
    end
    avrDE /= halfperiod
    avrDL /= halfperiod
    avrDEE /= halfperiod
    avrDEL /= halfperiod
    avrDLL /= halfperiod

    return avrDE, avrDL, avrDEE, avrDEL, avrDLL
end
