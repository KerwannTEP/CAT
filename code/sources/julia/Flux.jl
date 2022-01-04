##################################################
# Test star DF
##################################################

"""
    _HdHTest(a,b,c,d,x)

Computes the `H(a,b,c,d;x)` function used in the distribution function and its derivative.
Returns a tuple `(H, dH/dx, d2H/dx2)`.


# Arguments
- `x::Float64`: `x` parameter.
"""
function _HdHTest(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64)
    if (x <= 1)
        pref = 1/(GAMMA_ca*GAMMA_ad)
        HG   = _₂F₁(a+b,1+a-c,a+d,x)#tabHGInt[1](x) #
        HGp  = _₂F₁(a+b+1,a-c+2,a+d+1,x)#tabHGInt[3](x) #
        HGpp = _₂F₁(a+b+2,a-c+3,a+d+2,x)#tabHGInt[5](x) #
        H = pref*x^a*HG
        dH = pref*x^(a-1)*(a*HG + (a+b)*(1+a-c)*x*HGp/(a+d))
        d2H = (pref*x^(a-2)*(a*(a-1)*HG + (2*a)*(a+b)*(1+a-c)*x*HGp/(a+d)
               + ((a+b)*(1+a-c)/(a+d))*((a+b+1)*(2+a-c)/(a+d+1)) *x^2*HGpp))
        return H, dH, d2H
    else
        pref = 1/(GAMMA_db*GAMMA_bc)
        HG   = _₂F₁(a+b,1+b-d,b+c,1/x)#tabHGInt[2](1/x) #
        HGp  = _₂F₁(a+b+1,b-d+2,b+c+1,1/x)#tabHGInt[4](1/x) #
        HGpp = _₂F₁(a+b+2,b-d+3,b+c+2,1/x)#tabHGInt[6](1/x) #
        H = pref*x^(-b)*HG
        dH = pref*x^(-b-1)*((-b)*HG - (1/x)*HGp*(a+b)*(1+b-d)/(b+c))
        d2H = (pref*x^(-b-2)*(b*(b+1)*HG + (2*b+2)*(a+b)*(1+b-d)*(1/x)*HGp/(b+c)
               + ((a+b)*(1+b-d)/(b+c))*((a+b+1)*(2+b-d)/(b+c+1)) *(1/x^2)*HGpp))
        return H, dH, d2H
    end
end




"""
    _tFdF(tE,tL,[q=qCalc])

Gradients of the non-dimensional distribution function of a Plummer sphere.
Returns a tuple `(dtF/dtE, dtF/dtL, d2tF/dtE2, d2tF/dtEtL, d2tF/dtL2)`.

# Arguments
- `tE::Float64`: Reduced energy `tE = E/E0`.
- `tL::Float64`: Reduced energy `tL = L/L0`.
- `q ::Float64`: Anisotropy parameter.
"""
function _tdFTest(tE::Float64, tL::Float64, q::Float64=qCalc)
    if (tE <= 0.0 || tL <= 0.0) # If E or L are negative, the DF vanishes
        dE = 0.0
        dL = 0.0
        dEE = 0.0
        dEL = 0.0
        dLL = 0.0
        return dE, dL, dEE, dEL, dLL
    else
        x = tL^2/(2*tE)
        if (q == 0.0) # Isotropic case
            dE = 3/(PI^3) * (2*tE)^(5/2)
            dL = 0.0
            dEE = 15/(PI^3) * (2*tE)^(3/2)
            dEL = 0.0
            dLL = 0.0
            return dE, dL, dEE, dEL, dLL
        else
            if (q == 2.0)
                if (x <= 1.0)
                    pref = 18/(2*PI)^3
                    dE = pref * (2*tE-tL^2)^(1/2)
                    dL = -pref*tL * (2*tE-tL^2)^(1/2)
                    dEE = pref* (2*tE-tL^2)^(-1/2)
                    dEL = -pref*tL * (2*tE-tL^2)^(-1/2)
                    dLL = -pref*((2*tE-tL^2)^(1/2) - tL^2*(2*tE-tL^2)^(-1/2))
                    # here L ???
                    return dE, dL, dEE, dEL, dLL
                else
                    dE = 0.0
                    dL = 0.0
                    dEE = 0.0
                    dEL = 0.0
                    dLL = 0.0
                    return dE, dL, dEE, dEL, dLL
                end
            else
                pref = 3*GAMMA_6q/(2*(2*PI)^(5/2))
                pref1 = pref*tE^(5/2-q)
                pref2 = pref*tE^(3/2-q)
                H, dH, d2H = _HdHTest(0.0,q/2,9/2-q,1.0,x)
                dE = pref1*((7/2-q)*H - x*dH)
                dL = pref1*(tL*dH)
                dEE = pref2*((7/2-q)*(5/2-q)*H - (5-2*q)*x*dH + x^2*d2H)
                dEL = pref2*tL*((5/2-q)*dH - x*d2H)
                dLL = pref2*(tE*dH +tL^2*d2H)
                return dE, dL, dEE, dEL, dLL
            end
        end
    end
end


"""
    _dTestF(E,L,[q=qCalc])

Gradients of the distribution function of a Plummer sphere.
Normalized to the total mass of the cluster `_M`.
Returns a tuple `(dF/dE, dF/dL, d2F/dE2, d2F/dEL, d2F/dL2)`.

# Arguments
- `E::Float64`: Energy.
- `L::Float64`: Angular momentum.
- `q ::Float64`: Anisotropy parameter.
"""
function _dFTest(E::Float64, L::Float64, q::Float64=qCalc)
    tE = _tE(E)
    tL = _tL(L)
    tdE, tdL, tdEE, tdEL, tdLL = _tdFTest(tE,tL,q)
    dE = (_M*_F0/_E0)*tdE
    dL = (_M*_F0/_L0)*tdL
    dEE = (_M*_F0/(_E0^2))*tdEE
    dEL = (_M*_F0/(_E0*_L0))*tdEL
    dLL = (_M*_F0/(_L0^2))*tdLL
    return dE, dL, dEE, dEL, dLL
end



##################################################
# Compute numerical-gradient flux and its divergence
##################################################

"""
    divflux_NR_num(Jr,L,m_field,[eps=10^(-5)*_L0,nbK=nbK_default,nbAvr=nbAvr_default,m_test=m_field])

Computes the divergence of the flux in action space.

# Remarks
- `dF/dt = -div(flux)`.

# Arguments
- `Jr::Float64`: Radial action parameter.
- `L::Float64`: Angular momentum parameter.
- `m_field::Float64`: Mass of field star.
- `eps::Float64`: Finite difference used in numerical gradient.
- `nbK::Int64`: Sampling number of Rosenbluth integrals.
- `nbAvr::Int64`: Sampling number of orbit average.
- `m_test::Float64`: Mass of test star.
"""
function divflux_NR_num(Jr::Float64, L::Float64, m_field::Float64, eps::Float64=10^(-5)*_L0 ,
                    nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
                    m_test::Float64=m_field)


    Jr_p = Jr+eps
    Jr_m = Jr-eps

    L_p = L+eps
    L_m = L-eps

    E = _E_from_Jr(Jr,L)
    E_Lp = _E_from_Jr(Jr,L_p,nbAvr)
    E_Lm = _E_from_Jr(Jr,L_m,nbAvr)
    E_p = _E_from_Jr(Jr_p,L,nbAvr)
    E_m = _E_from_Jr(Jr_m,L,nbAvr)
    E_p_Lp = _E_from_Jr(Jr_p,L_p,nbAvr)
    E_p_Lm = _E_from_Jr(Jr_p,L_m,nbAvr)
    E_m_Lp = _E_from_Jr(Jr_m,L_p,nbAvr)
    E_m_Lm = _E_from_Jr(Jr_m,L_m,nbAvr)

    # Distribution function
    Ftot = _F(E,L)
    dFtotdE, dFtotdL, d2FtotdE2, d2FtotdEL, d2FtotdL2 = _dFTest(E,L)

    # Compute dJr/dE and dJr/dL
    sp, sa = radius_s_bounds_action(Jr,L,nbAvr)
    a = sma_eff(sp,sa)
    ecc = ecc_eff(sp,sa)

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

        dJrdE += jac_loc
        dJrdL += jac_loc/rloc^2
        d2JrdE2 += djacdE_loc
        d2JrdEL += djacdL_loc
        d2JrdL2 += (jac_loc/trloc^2+L/trloc^4*
                (djacdL_loc*trloc^2-2*jac_loc*sloc*dsdL_loc))

    end

    dJrdE *= (2/nbAvr)*(1/_Omega0)/PI
    dJrdL *= (2/nbAvr)*(-L/PI)*(1/_Omega0)
    d2JrdE2 *= (2/nbAvr)*(1/_Omega0)/PI
    d2JrdEL *= (2/nbAvr)*(1/_Omega0)/PI
    d2JrdL2 *= (2/nbAvr)*(-1/PI)/_L0

    # Compute DFtot/DJr and DFtot/DL
    dEdL = -dJrdL/dJrdE
    dEdJr = 1/dJrdE

    d2EdJr2 = -d2JrdE2/(dJrdE)^3
    d2EdJrL = -(d2JrdE2*dEdL+d2JrdEL)/(dJrdE)^2
    d2EdL2 = -((dEdL*d2JrdE2+d2JrdEL)*dEdL+dEdL*d2JrdEL+d2JrdL2)/(dJrdE)

    DFtotDJr = (1/dJrdE)*dFtotdE
    DFtotDL = dEdL*dFtotdE + dFtotdL

    D2FtotDJr2 = d2EdJr2*dFtotdE+(dEdJr)^2*d2FtotdE2
    D2FtotDJrL = d2EdJrL*dFtotdE + dEdJr*(dEdL*d2FtotdE2+d2FtotdEL)
    D2FtotDL2 = d2EdL2*dFtotdE + dEdL*(dEdL*d2FtotdE2+d2FtotdEL)+dEdL*d2FtotdEL+d2FtotdL2

    # Diffusion coefficients
    DJr, DL, DJrJr, DJrL, DLL = avr_action_coefficients(Jr,L,m_field,nbK,nbAvr,m_test)
    DJr_p, _, DJrJr_p, DJrL_Jrp, _ = avr_action_coefficients(Jr_p,L,m_field,nbK,nbAvr,m_test)
    DJr_m, _, DJrJr_m, DJrL_Jrm, _ = avr_action_coefficients(Jr_m,L,m_field,nbK,nbAvr,m_test)
    _, DL_p, _, DJrL_Lp, DLL_p = avr_action_coefficients(Jr,L_p,m_field,nbK,nbAvr,m_test)
    _, DL_m, _, DJrL_Lm, DLL_m = avr_action_coefficients(Jr,L_m,m_field,nbK,nbAvr,m_test)
    _, _, _, DJrL_pp, _ = avr_action_coefficients(Jr_p,L_p,m_field,nbK,nbAvr,m_test)
    _, _, _, DJrL_mm, _ = avr_action_coefficients(Jr_m,L_m,m_field,nbK,nbAvr,m_test)
    _, _, _, DJrL_pm, _ = avr_action_coefficients(Jr_p,L_m,m_field,nbK,nbAvr,m_test)
    _, _, _, DJrL_mp, _ = avr_action_coefficients(Jr_m,L_p,m_field,nbK,nbAvr,m_test)


    # d(DJrF)/dJr
    termJr1 = (2*L*DJr_p-2*L*DJr_m)/(2*eps*_L0)*Ftot
    termJr2 = 2*L*DJr*DFtotDJr
    termJr = termJr1 + termJr2

    # d(DLF)/dL
    termL1 = (2*L_p*DL_p-2*L_m*DL_m)/(2*eps)*Ftot
    termL2 = 2*L*DL*DFtotDL
    termL = termL1 + termL2

    # d2(DJrJrF)/dJr2
    termJrJr1 = (2*L*DJrJr_p+2*L*DJrJr_m-2*2*L*DJrJr)/((eps)^2)*Ftot
    termJrJr2 = 2*(2*L*DJrJr_p-2*L*DJrJr_m)/(2*eps)*DFtotDJr
    termJrJr3 = 2*L*DJrJr*D2FtotDJr2
    termJrJr = termJrJr1 + termJrJr2 + termJrJr3

    # d2(DLLF)/dL2
    termLL1 = (2*L_p*DLL_p+2*L_m*DLL_m-2*2*L*DLL)/((eps)^2)*Ftot
    termLL2 = 2*(2*L_p*DLL_p-2*L_m*DLL_m)/(2*eps)*DFtotDL
    termLL3 = 2*L*DLL*D2FtotDL2
    termLL = termLL1 + termLL2 + termLL3

    # d2(DJrLF)/dJrL
    termJrL1 = (2*L_p*DJrL_pp+2*L_m*DJrL_mm-2*L_m*DJrL_pm-2*L_p*DJrL_mp)/(4*(eps)^2)*Ftot
    termJrL2 = (2*L*DJrL_Jrp-2*L*DJrL_Jrm)/(2*eps)*DFtotDL
    termJrL3 = (2*L_p*DJrL_Lp-2*L_m*DJrL_Lm)/(2*eps)*DFtotDJr
    termJrL4 = 2*L*DJrL*D2FtotDJrL
    termJrL = termJrL1 + termJrL2 + termJrL3 + termJrL4

    return termJr+termL-0.5*termJrJr-0.5*termLL-termJrL
end

"""
    flux_NR_num(Jr,L,m_field,[eps=10^(-5)*_L0,nbK=nbK_default,nbAvr=nbAvr_default,m_test=m_field])

Computes the flux in action space.

# Arguments
- `Jr::Float64`: Radial action parameter.
- `L::Float64`: Angular momentum parameter.
- `m_field::Float64`: Mass of field star.
- `eps::Float64`: Finite difference used in numerical gradient.
- `nbK::Int64`: Sampling number of Rosenbluth integrals.
- `nbAvr::Int64`: Sampling number of orbit average.
- `m_test::Float64`: Mass of test star.
"""
function flux_NR_num(Jr::Float64, L::Float64, m_field::Float64, eps::Float64=10^(-5)*_L0,
                    nbK::Int64=nbK_default, nbAvr::Int64=nbAvr_default,
                    m_test::Float64=m_field)

        Jr_p = Jr+eps
        Jr_m = Jr-eps

        L_p = L+eps
        L_m = L-eps

        E = _E_from_Jr(Jr,L)

        Ftot, dFtotdE, dFtotdL = _FdF(E,L)

        # Compute dJr/dE and dJr/dL
        sp, sa = radius_s_bounds_action(Jr,L,nbAvr)
        a = sma_eff(sp,sa)
        ecc = ecc_eff(sp,sa)

        dJrdE = 0.0
        dJrdL = 0.0

        for iu=1:nbAvr
            uloc = -1+2*(iu-0.5)/nbAvr
            rloc = r_eff_anomaly(uloc,a,ecc)

            jac_loc = dth1du(uloc,sp,sa)
            dJrdE += jac_loc
            dJrdL += jac_loc/rloc^2
        end

        dJrdE *= (2/nbAvr)*(1/_Omega0)/PI
        dJrdL *= (2/nbAvr)*(-L/PI)*(1/_Omega0)

        # Compute DFtot/DJr and DFtot/DL
        DFtotDJr = (1/dJrdE)*dFtotdE

        dEdL = -dJrdL/dJrdE
        DFtotDL = dEdL*dFtotdE + dFtotdL

        # Compute diffusion coefficients
        DJr, DL, DJrJr, DJrL, DLL = avr_action_coefficients(Jr,L,m_field,nbK,nbAvr,m_test) #
        _, _, DJrJr_Jrp, DJrL_Jrp, _ = avr_action_coefficients(Jr_p,L,m_field,nbK,nbAvr,m_test) #
        _, _, DJrJr_Jrm, DJrL_Jrm, _ = avr_action_coefficients(Jr_m,L,m_field,nbK,nbAvr,m_test) #
        _, _, _, DJrL_Lp, DLL_Lp = avr_action_coefficients(Jr,L_p,m_field,nbK,nbAvr,m_test) #
        _, _, _, DJrL_Lm, DLL_Lm = avr_action_coefficients(Jr,L_m,m_field,nbK,nbAvr,m_test) #

        # Flux Jr
        termJr1 = 2*L*DJr*Ftot
        termJr2 = -0.5*((2*L*DJrJr_Jrp-2*L*DJrJr_Jrm)/(2*eps))*Ftot
        termJr3 = -L*DJrJr*DFtotDJr
        termJr4 = -0.5*((2*L_p*DJrL_Lp-2*L_m*DJrL_Lm)/(2*eps))*Ftot
        termJr5 = -L*DJrL*DFtotDL

        FJr = termJr1+termJr2+termJr3+termJr4+termJr5

        # Flux Jr
        termL1 = 2*L*DL*Ftot
        termL2 = -0.5*((2*L*DJrL_Jrp-2*L*DJrL_Jrm)/(2*eps))*Ftot
        termL3 = -L*DJrL*DFtotDJr
        termL4 = -0.5*((2*L_p*DLL_Lp-2*L_m*DLL_Lm)/(2*eps))*Ftot
        termL5 = -L*DLL*DFtotDL

        FL = termL1+termL2+termL3+termL4+termL5

        return FJr, FL

end
