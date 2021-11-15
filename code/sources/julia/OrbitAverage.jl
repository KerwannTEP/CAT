##################################################
# Orbit-averagering of the diffusion coefficients
# In (E,L) space
##################################################
using PolynomialRoots

"""
    radius_s_bounds(E,L)

Computes the bounds of the orbit with action parameters `(E,L)`.
Bounds are in the variable `s`, defined as `s^2 = 1 + x^2`, where `x = r/b`.

Returns `"Unbounded orbit"` for unbounded orbits.
Returns `"Not a possible orbit"` for forbidden parameters.

# Remarks
- Uses `PolynomialRoots`.

# Arguments
- `E::Float64`: Energy parameter.
- `L ::Float64`: Angular momentum parameter.
"""
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

"""
    sma_eff(sp,sa)

Computes the effective semi-major axis of an orbit with `s`-bounds `sp` and `sa`.
Bounds are in the variable `s`, defined as `s^2 = 1 + x^2`, where `x = r/b`.

# Arguments
- `sp::Float64`: Pericenter in the radius variable `s`.
- `sa::Float64`: Apocenter in the radius variable `s`.
"""
function sma_eff(sp::Float64, sa::Float64)
    return (sa+sp)/2
end

"""
    ecc_eff(sp,sa)

Computes the effective eccentricity of an orbit with `s`-bounds `sp` and `sa`.
Bounds are in the variable `s`, defined as `s^2 = 1 + x^2`, where `x = r/b`.

# Arguments
- `sp::Float64`: Pericenter in the radius variable `s`.
- `sa::Float64`: Apocenter in the radius variable `s`.
"""
function ecc_eff(sp::Float64, sa::Float64)
    return (sa-sp)/(sa+sp)
end

"""
    tr_eff_anomaly(u,a,ecc)

Computes the non-dimensional radius of a star with effective anomaly `u`, semi-major axis `a` and eccentricity `ecc`.

# Arguments
- `u::Float64`: Effective anomaly parameter.
- `a::Float64`: Effective semi-major axis.
- `ecc::Float64`: Effective eccentricity.
"""
function tr_eff_anomaly(u::Float64, a::Float64, ecc::Float64)
    s = a*(1+ecc*u*(1.5-u^2/2))

    return sqrt(abs(s^2-1)) # use absolute value for when s->1
end

"""
    r_eff_anomaly(u,a,ecc)

Computes the radius of a star with effective anomaly `u`, semi-major axis `a` and eccentricity `ecc`.

# Arguments
- `u::Float64`: Effective anomaly parameter.
- `a::Float64`: Effective semi-major axis.
- `ecc::Float64`: Effective eccentricity.
"""
function r_eff_anomaly(u::Float64, a::Float64, ecc::Float64)
    return _b*tr_eff_anomaly(u,a,ecc)
end

"""
    dth1du(u,sp,sa)

Non-dimensional Jacobian of the transformation from angles `θ` to effective anomaly `u`.

# Arguments
- `u::Float64`: Effective anomaly.
- `sp::Float64`: Pericenter in the variable `s`.
- `sa::Float64`: Apocenter in the variable `s`.
"""
function dth1du(u::Float64, sp::Float64, sa::Float64)
    num = sa*(-2+u)*(1+u)^2-sp*(2-3*u+u^3)
    denumSq = ((sa^2*sp*(-2+u)*(1+u)^2-sp*(6-3*u+u^3)
               +sa*(-6-3*u+u^3-sp^2*(2-3*u+u^3)))/
               (sa*sp*(sa+sp)*(sa*(-2+u)*(1+u)^2-sp*(2-3*u+u^3))))
    return -3/(4*sqrt(2))*1/sqrt(4-u^2)*num/sqrt(denumSq)
end

"""
    halfperiod(E,L,[nbAvr=nbAvr_default])

Computes half the radial period of the orbit with parameters `(E,L)`.

# Remarks
- Well-defined everywhere by using an effective anomaly.

# Arguments
- `E::Float64`: Energy parameter.
- `L::Float64`: Angular momentum parameter.
- `nbAvr::Int64`: Orbit average sampling number.
"""
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

"""
    avr_action_coefficients(E,L,m_field,[nbK=nbK_default,nbAvr=nbAvr_default,m_test=m_field])

Orbit-averaged diffusion coefficients in energy/angular momentum space `(E,L)`.
Returns a 5-uple `(avrDE, avrDL, avrDEE, avrDEL, avrDLL)`.

# Arguments
- `E::Float64`: Energy parameter.
- `L::Float64`: Angular momentum parameter.
- `m_field::Float64`: Mass of field star.
- `nbK::Int64`: Sampling number of Rosenbluth integrals.
- `nbAvr::Int64`: Sampling number of orbit average.
- `m_test::Float64`: Mass of test star.
"""
function avr_orbit_coefficients(E::Float64, L::Float64, m_field::Float64,
                                nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
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
