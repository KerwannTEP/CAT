##################################################
# Go into (Jr,L) space
##################################################

# more adapted version
# always takes possible orbital parameters
function radius_s_bounds_action(Jr::Float64, L::Float64, nbAvr::Int64=nbAvr_default)

    E = _E_from_Jr(Jr,L,nbAvr)

    tE = _tE(E)
    tL = _tL(L)
    if (L != 0.0)
        rts = sort(real(roots([1,-(tE-tL^2/2),-1,tE])))
        return Float64(rts[2]), Float64(rts[3])
    else
        return 1.0, 1/tE
    end

end

function halfperiod_action(Jr::Float64, L::Float64, nbAvr::Int64=nbAvr_default)

    sp, sa = radius_s_bounds_action(Jr,L)
    a = sma_eff(sp,sa)
    ecc = ecc_eff(sp,sa)

    halfperiod = 0.0

    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        rloc = r_eff_anomaly(uloc,a,ecc)
        jac_loc = dth1du(uloc,sp,sa)

        halfperiod += jac_loc
    end
    halfperiod *= (1/_Omega0)*(2/nbAvr)

    return halfperiod
end


function PHIeff(s::Float64, L::Float64)
    return _E0/s + L^2/(2*_b^2*(s^2-1))
end

# returns dPHIeff/ds and dPHIeff/dL
function dPHIeff(s::Float64, L::Float64)
    ds = -_E0/s^2 - s*L^2/(_b^2*(s^2-1)^2)
    dL = L/(_b^2*(s^2-1))
    return ds, dL
end

function _fAnomaly(u::Float64)
    return u*(1.5-u^2/2)
end

function _djacdsp(u::Float64, sp::Float64, sa::Float64)

    num = (3 *(4 *sp^3 *(-1 + u)^2 *(12 - 3 *u^2 + 2 *u^3 + u^4) +
    sa *sp^2 *(108 - 120 *u - 27 *u^2 + 40 *u^3 + 18 *u^4 - 3 *u^6 +
       3 *sp^2 *(-1 + u)^4 *(2 + u)^2) -
    2 *sa^2 *sp *(-6 - 3 *u + u^3) *(6 - 3 *u + u^3 +
       sp^2 *(2 - 3 *u + u^3)) -
    sa^3 *(-2 + u) *(1 + u)^2 *(6 + 3* u - u^3 + sp^2 *(6 - 3 *u + u^3))))

    den = (8 *sqrt(2)*
   sp *(sa +
    sp) *sqrt(-(((-4 + u^2) *(sa^2 *sp *(-2 + u) *(1 + u)^2 -
      sp *(6 - 3 *u + u^3) +
      sa *(-6 - 3 *u + u^3 - sp^2 *(2 - 3* u + u^3))))/(
   sa *sp *(sa + sp) *(sa *(-2 + u) *(1 + u)^2 -
      sp *(2 - 3 *u + u^3)))))* (sa^2* sp *(2 + 3 *u - u^3) +
    sp *(6 - 3* u + u^3) + sa *(6 + 3*u - u^3 + sp^2 *(2 - 3 *u + u^3))))

    return num/den
end

function _djacdsa(u::Float64, sp::Float64, sa::Float64)

    num = (3 *(-4 + u^2) *(3 *sa^4 *sp *(-2 + u)^2 *(1 + u)^4 +
    sp^3 *(-1 + u)^2 *(12 - 3 *u^2 + 2 *u^3 + u^4) -
    2 *sa *sp^2 *(-36 + 9 *u^2 - 6 *u^4 + u^6) -
    2 *sa^3 *(-2 + u) *(1 + u)^2 *(12 + 6 *u - 2 *u^3 +
       sp^2 *(6 - 3 *u + u^3)) +
    sa^2 *sp *(108 + 120 *u - 27 *u^2 - 40*u^3 + 18 *u^4 - 3 *u^6 -
       sp^2 *(-12 + 12 *u + 9 *u^2 - 4 *u^3 - 6 *u^4 + u^6))))

    den =(8 *sqrt(2)*
   sa^2 *sp *(sa + sp)^2 *(sa *(-2 + u) *(1 + u)^2 -
    sp *(2 - 3 *u + u^3)) *(-(((-4 + u^2) *(sa^2 *sp *(-2 + u) *(1 + u)^2 -
       sp *(6 - 3* u + u^3) +
       sa *(-6 - 3* u + u^3 - sp^2 *(2 - 3 *u + u^3))))/(
    sa *sp* (sa + sp) *(sa *(-2 + u) *(1 + u)^2 - sp *(2 - 3 *u + u^3)))))^(
  3/2))

    return num/den
end

# returns djacdE, djacdL, dsdL
function djac_and_ds(u::Float64, sp::Float64, sa::Float64, L::Float64)
    dPHIdsp, dPHIdLp = dPHIeff(sp,L)
    dPHIdsa, dPHIdLa = dPHIeff(sa,L)
    dspdE = 1/dPHIdsp
    dsadE = 1/dPHIdsa
    dspdL = -dPHIdLp/dPHIdsp
    dsadL = -dPHIdLa/dPHIdsa

    djacdsp = _djacdsp(u,sp,sa)
    djacdsa = _djacdsa(u,sp,sa)

    dsdsp = (1-_fAnomaly(u))/2
    dsdsa = (1+_fAnomaly(u))/2


    djacdE = dsadE*djacdsa + dspdE*djacdsp
    djacdL = dsadL*djacdsa + dspdL*djacdsp

    dsdL = dsadL*dsdsa+dspdL*dsdsp

    return djacdE, djacdL, dsdL
end

function vrSq(r::Float64, E::Float64, L::Float64)
    return 2*(E-psiEff(r,L))
end


#function avr_action_coefficients(E::Float64, L::Float64, m_field::Float64,
function avr_action_coefficients(Jr::Float64, L::Float64, m_field::Float64,
                                nbK::Int=nbK_default, nbAvr::Int64=nbAvr_default,
                                m_test::Float64=m_field)

    E = _E_from_Jr(Jr,L,nbAvr)
    sp, sa = radius_s_bounds_action(Jr,L,nbAvr) #radius_s_bounds(E,L)


    a = sma_eff(sp,sa)
    ecc = ecc_eff(sp,sa)

    avrDE = 0.0
    avrDL = 0.0
    avrDEE = 0.0
    avrDEL = 0.0
    avrDLL = 0.0
    halfperiod = 0.0

    dJrdE = 0.0
    dJrdL = 0.0
    d2JrdE2 = 0.0
    d2JrdEL = 0.0
    d2JrdL2 = 0.0


    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        rloc = r_eff_anomaly(uloc,a,ecc)
        trloc = rloc/_b
        sloc = sqrt(trloc^2+1)

        jac_loc = dth1du(uloc,sp,sa)
        djacdE_loc, djacdL_loc, dsdL_loc = djac_and_ds(uloc,sp,sa,L)

        dEloc, dE2loc, dLloc, dEdLloc, dL2loc = localOrbitChange(rloc,E,L,m_field,nbK,m_test)
        vrSqloc = vrSq(rloc,E,L)

        avrDE += jac_loc*dEloc
        avrDL += jac_loc*dLloc
        avrDEE += jac_loc*dE2loc
        avrDEL += jac_loc*dEdLloc
        avrDLL += jac_loc*dL2loc
        halfperiod += jac_loc

        dJrdE += jac_loc
        dJrdL += jac_loc/rloc^2
        d2JrdE2 += djacdE_loc
        d2JrdEL += djacdL_loc
        d2JrdL2 += (jac_loc/trloc^2+L/trloc^4*
                (djacdL_loc*trloc^2-2*jac_loc*sloc*dsdL_loc))

    end

    avrDE /= halfperiod
    avrDL /= halfperiod
    avrDEE /= halfperiod
    avrDEL /= halfperiod
    avrDLL /= halfperiod

    dJrdE *= (2/nbAvr)*(1/_Omega0)/PI
    dJrdL *= (2/nbAvr)*(-L/PI)*(1/_Omega0)
    d2JrdE2 *= (2/nbAvr)*(1/_Omega0)/PI
    d2JrdEL *= (2/nbAvr)*(1/_Omega0)/PI
    d2JrdL2 *= (2/nbAvr)*(-1/PI)/_L0

    avrDJr = dJrdE*avrDE+dJrdL*avrDL+(1/2)*d2JrdE2*avrDEE+(1/2)*d2JrdL2*avrDLL+d2JrdEL*avrDEL
    avrDJrJr = dJrdE^2*avrDEE+dJrdL^2*avrDLL+2*dJrdE*dJrdL*avrDEL
    avrDJrL = dJrdE*avrDEL+dJrdL*avrDLL

    return avrDJr, avrDL, avrDJrJr, avrDJrL, avrDLL
end


function _Jr(E::Float64, L::Float64, nbAvr::Int64=100)
    Jr = 0.0

    sp, sa = radius_s_bounds(E,L)
    a = sma_eff(sp,sa)
    ecc = ecc_eff(sp,sa)

    for iu=1:nbAvr
        uloc = -1+2*(iu-0.5)/nbAvr
        rloc = r_eff_anomaly(uloc,a,ecc)
        jac_loc = dth1du(uloc,sp,sa)
        vrSqloc = vrSq(rloc,E,L)

        Jr += jac_loc*vrSqloc
    end

    Jr *= (2/nbAvr)*(1/_Omega0)/PI
    return Jr
end

# returns E(Jr,L)
# unstable for very small L>0
function _E_from_Jr(Jr::Float64, L::Float64, nbAvr::Int64=nbAvr_default, prec::Float64=4.0*10^(-15))
    En = -0.5 # E0
    delta = 0.5 # initial precision
    JrE = 0.0
    if (Jr==0.0 && L==0.0)
        return -1.0
    else
        while (delta>prec) # while E is not a solution precision enough
            delta /= 2
            if (L > Lc(_tE(En))) # if (En,L) is not an orbit, we move towards the accepted (E,L) region
                En += delta
            else # if within the accepted (E,L) region, use the fact that Jr(E,L) increases with E
                JrE = _Jr(En,L,nbAvr)
                if (JrE < Jr)
                    En += delta
                else
                    En -= delta
                end
            end
        end
        return En
    end
end


function reduced_DF(Jr::Float64, L::Float64,nbAvr::Int64=nbAvr_default)
    E = _E_from_Jr(Jr,L,nbAvr)
    return 2*L*_F(E,L)
end
