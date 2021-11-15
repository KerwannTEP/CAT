##################################################
# Useful general functions
##################################################
using StaticArrays # To have access to static arrays

using SpecialFunctions
using HypergeometricFunctions

##################################################
# Compute normalized variables
##################################################

"""
    _tr(r)

Computes the non-dimensional radius `x = r/b`.

# Arguments
- `r::Float64`: Radius.
"""
function _tr(r::Float64)
    return r/_b
end

"""
    _tE(E)

Computes the non-dimensional energy `tE = E/E0`.

# Arguments
- `E::Float64`: Energy.
"""
function _tE(E::Float64)
    return E/_E0
end

"""
    _tL(L)

Computes the non-dimensional angular momentum `tL = L/L0`.

# Arguments
- `L::Float64`: Angular momentum.
"""
function _tL(L::Float64)
    return L/_L0
end

##################################################
# Compute the normalized Plummer potential
##################################################

"""
    _tpsi(tr)

Computes the non-dimensional Plummer potential.

# Arguments
- `tr::Float64`: Non-dimensional radius.
"""
function _tpsi(tr::Float64)
    return 1/sqrt(1.0+tr^2)
end

"""
    _tdpsidr(tr)

Computes the derivative of the non-dimensional Plummer potential.

# Arguments
- `tr::Float64`: Non-dimensional radius.
"""
function _tdpsidr(tr::Float64)
    return -tr/(1.0+tr^2)^(3/2)
end

"""
    _td2psidr2(tr)

Computes the second derivative of the non-dimensional Plummer potential.

# Arguments
- `tr::Float64`: Non-dimensional radius.
"""
function _td2psidr2(tr::Float64)
    return -(1.0-tr^2)/(1.0+tr^2)^(5/2)
end

##################################################
# Compute the normalized effective Plummer potential
##################################################

"""
    _tpsiEff(tr,tL)

Computes the non-dimensional effective Plummer potential.

# Arguments
- `tr::Float64`: Non-dimensional radius.
- `tL::Float64`: Non-dimensional angular momentum.
"""
function _tpsiEff(tr::Float64, tL::Float64)
    return _tpsi(tr) - tL^2/(2*tr^2)
end

"""
    _tdpsiEffdr(tr,tL)

Computes the tr-derivative of the non-dimensional effective Plummer potential.

# Arguments
- `tr::Float64`: Non-dimensional radius.
- `tL::Float64`: Non-dimensional angular momentum.
"""
function _tdpsiEffdr(tr::Float64, tL::Float64)
    return _tdpsidr(tr) + tL^2/(tr^3)
end

"""
    _td2psiEffdr2(tr,tL)

Computes the second tr-derivative of the non-dimensional effective Plummer potential.

# Arguments
- `tr::Float64`: Non-dimensional radius.
- `tL::Float64`: Non-dimensional angular momentum.
"""
function _td2psiEffdr2(tr::Float64, tL::Float64)
    return _td2psidr2(tr) - 3* tL^2/(tr^4)
end


##################################################
# Compute the Plummer potential
##################################################

"""
    psi(r)

Computes the Plummer potential.

# Arguments
- `r::Float64`: Radius.
"""
function psi(r::Float64)
    return _E0*_tpsi(_tr(r))
end

"""
    dpsidr(r)

Computes the derivative of the Plummer potential.

# Arguments
- `r::Float64`: Radius.
"""
function dpsidr(r::Float64)
    return (_E0/_b)*_tdpsidr(_tr(r))
end

"""
    d2psidr2(r)

Computes the second derivative of the Plummer potential.

# Arguments
- `r::Float64`: Radius.
"""
function d2psidr2(r::Float64)
    return (_E0/_b^2)*_td2psidr2(_tr(r))
end

##################################################
# Compute the effective Plummer potential
##################################################

"""
    psiEff(r,L)

Computes the effective Plummer potential.

# Arguments
- `r::Float64`: Radius.
- `L::Float64`: Angular momentum.
"""
function psiEff(r::Float64, L::Float64)
    return _E0*_tpsiEff(_tr(r),_tL(L))
end

"""
    dpsiEffdr(r,L)

Computes the r-derivative of the effective Plummer potential.

# Arguments
- `r::Float64`: Radius.
- `L::Float64`: Angular momentum.
"""
function dpsiEffdr(r::Float64, L::Float64)
    return (_E0/_b)*_tdpsiEffdr(_tr(r),_tL(L))
end

"""
    d2psiEffdr2(r,L)

Computes the second r-derivative of the effective Plummer potential.

# Arguments
- `r::Float64`: Radius.
- `L::Float64`: Angular momentum.
"""
function d2psiEffdr2(r::Float64, L::Float64)
    return (_E0/_b^2)*_td2psiEffdr2(_tr(r),_tL(L))
end

##################################################
# Compute the Plummer density
##################################################

"""
    _trho(tr)

Computes the non-dimensional Plummer density.

# Arguments
- `tr::Float64`: Non-dimensional radius.
"""
function _trho(tr::Float64)
    return 1/(1+tr^2)^(5/2)
end

"""
    rho(r)

Computes the Plummer density.

# Arguments
- `r::Float64`: Radius.
"""
function rho(r::Float64)
    return _rho0*_trho(_tr(r))
end

##################################################
# Compute Lc(E) the angular momentum per unit mass of a circular orbit
##################################################

"""
    _tEta(s,tE)

Function used to compute the maximal, circular angular momentum of an orbit with non-dimensional energy `tE`.
Uses the radius parameter `s` defined by `s^2 = 1 + x^2` where `x = r/b`.

# Arguments
- `s::Float64`: Radius parameter.
- `tE::Float64`: Non-dimensional energy.
"""
function _tEta(s::Float64, tE::Float64)
    return (s^2-1)*(1/s - tE)
end

"""
    _sc(tE)

Compute the s-radius of the circular orbit with non-dimensional energy `tE`.
Uses the radius parameter `s` defined by `s^2 = 1 + x^2` where `x = r/b`.

# Remarks
- Uses an analytical expression (cubic root).

# Arguments
- `tE::Float64`: Non-dimensional energy.
"""
function _sc(tE::Float64)
    t1 = 1/(6*tE)
    t2 = (1+54*tE^2)/(216*tE^3)
    t3 = (1/(4*tE))*sqrt(1+1/(27*tE^2))
    return t1+cbrt(t2+t3)+cbrt(t2-t3)
end

"""
    _rc(tE)

Compute the radius of the circular orbit with non-dimensional energy `tE`.

# Remarks
- Uses an analytical expression (cubic root).

# Arguments
- `tE::Float64`: Non-dimensional energy.
"""
function _rc(tE::Float64)
    sc = _sc(tE)
    return _b*sqrt(sc^2-1)
end

"""
    Lc(tE)

Compute the angular momentum of the circular orbit with non-dimensional energy `tE`.

# Remarks
- Uses an analytical expression (cubic root).

# Arguments
- `tE::Float64`: Non-dimensional energy.
"""
function Lc(tE::Float64)
    if (tE == 1.0)
        return 0.0
    else
        return _L0*sqrt(2*_tEta(_sc(tE),tE))
    end
end

##################################################
# Convert (vr,vt) to (E,L) at a given radius r
##################################################

"""
    energy(r,vr,vt)

Compute the energy of a star at position `r` and velocity `(vr,vt)`.

# Arguments
- `r::Float64`: Radius.
- `vr::Float64`: Radial velocity.
- `vt::Float64`: Tangential velocity.
"""
function energy(r::Float64, vr::Float64, vt::Float64)
    return psi(r) + vr^2/(2.0) + vt^2/(2.0)
end

"""
    momentum(r,vr,vt)

Compute the angular momentum of a star at position `r` and velocity `(vr,vt)`.

# Arguments
- `r::Float64`: Radius.
- `vr::Float64`: Radial velocity.
- `vt::Float64`: Tangential velocity.
"""
function momentum(r::Float64, vr::Float64, vt::Float64)
    return r*vt
end

##################################################
# Convert (E,L) to (vr,vt) at a given radius r
# Convention vr >= 0
##################################################

"""
    radialVelocitySq(r,E,L)

Compute the square radial velocity of a star at position `r` of an orbit with parameters `(E,L)`.

# Arguments
- `r::Float64`: Radius.
- `E::Float64`: Energy.
- `L::Float64`: Angular momentum.
"""
function radialVelocitySq(r::Float64,E::Float64, L::Float64)
    return 2*(E-psiEff(r,L))
end

"""
    radialVelocity(r,E,L)

Compute the radial velocity of a star at position `r` of an orbit with parameters `(E,L)`.
Uses the convention `vr >= 0`.

# Arguments
- `r::Float64`: Radius.
- `E::Float64`: Energy.
- `L::Float64`: Angular momentum.
"""
function radialVelocity(r::Float64,E::Float64, L::Float64)
    return sqrt(2*abs((E-psiEff(r,L))))
end

"""
    tangentVelocity(r,E,L)

Compute the tangential velocity of a star at position `r` of an orbit with parameters `(E,L)`.

# Arguments
- `r::Float64`: Radius.
- `E::Float64`: Energy.
- `L::Float64`: Angular momentum.
"""
function tangentVelocity(r::Float64,E::Float64, L::Float64)
    return L/r
end

"""
    velocity(r,E,L)

Compute the velocity of a star at position `r` of an orbit with parameters `(E,L)`.

# Arguments
- `r::Float64`: Radius.
- `E::Float64`: Energy.
- `L::Float64`: Angular momentum.
"""
function velocity(r::Float64,E::Float64, L::Float64)
    return sqrt(2*abs((E-psi(r,L))))
end



##################################################
# Bisection algorithm to find zeroes of functions
# Here, used to compute Ec(L)
##################################################

"""
    bisection(fun, xl, xu, <keyword arguments>)

Algorithm which looks for zeroes of the functions `fun` given in argument.
It is a simple bisection algorithm, but it makes no allocations and is sufficiently fast.
It allows us not to have to use the Roots library that makes a lot of allocations.
The optional tolerances are set to the same as the ones found in Roots.jl, and are most likely overkill.

The zero is searched within the interval `[xl, xu]`, whose bounds are given as arguments.
"""
function bisection(fun, xl::Float64, xu::Float64, tolx::Float64=4.0*10^(-15), tolf::Float64=4.0*eps(Float64), iterMAX::Int64=50)
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####
    if (abs(fl) <= tolf) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (abs(fu) <= tolf) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end
    #####
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET"
    #####
    iter = 0 # Counter for the iterations
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if ((abs(xu-xl) <= tolx) || (iter > iterMAX)) # The considered bracket is smaller than the tolerance, or we have made too many iterations
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        iter += 1 # Updating the counter of iterations
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
end


"""
    _rc_L(L)

Compute the radius of the circular orbit with angular momentum `L`.

# Remarks
- Uses a bisection algorithm.

# Arguments
- `L::Float64`: Angular momentum.
"""
function _rc_L(L::Float64)
    tL = _tL(L)
    fct = (tr->_tdpsiEffdr(tr,tL))
    trmin = tL^(2/3) # derivative is positive
    trmax = trmin
    while (_tdpsiEffdr(trmax,tL) > 0) # find value of trmin such that derivative is negative
        trmax *= 2
    end

    troot = bisection(fct,trmin,trmax)

    return _b*troot
end

"""
    Ec(L)

Compute the energy of the circular orbit with angular momentum `L`.

# Remarks
- Uses a bisection algorithm.

# Arguments
- `L::Float64`: Angular momentum.
"""
function Ec(L::Float64)
    rc = _rc_L(L)
    return psiEff(rc,L)
end
