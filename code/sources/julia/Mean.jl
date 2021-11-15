##################################################
# Useful general functions
##################################################
using StaticArrays # To have access to static arrays

using SpecialFunctions
using HypergeometricFunctions

##################################################
# Compute normalized quantites
##################################################

function _tr(r::Float64)
    return r/_b
end

function _tE(E::Float64)
    return E/_E0
end

function _tL(L::Float64)
    return L/_L0
end

function _tpsi(tr::Float64)
    return 1/sqrt(1.0+tr^2)
end

function _tdpsidr(tr::Float64)
    return -tr/(1.0+tr^2)^(3/2)
end

function _td2psidr2(tr::Float64)
    return -(1.0-tr^2)/(1.0+tr^2)^(5/2)
end

function _tpsiEff(tr::Float64, tL::Float64)
    return _tpsi(tr) - tL^2/(2*tr^2)
end

function _tdpsiEffdr(tr::Float64, tL::Float64)
    return _tdpsidr(tr) + tL^2/(tr^3)
end

function _td2psiEffdr2(tr::Float64, tL::Float64)
    return _td2psidr2(tr) - 3* tL^2/(tr^4)
end

##################################################
# Compute Ec(L) the binding energy per unit mass of a circular orbit
##################################################

# Potential
function psi(r::Float64)
    return _E0*_tpsi(_tr(r))
end

function dpsidr(r::Float64)
    return (_E0/_b)*_tdpsidr(_tr(r))
end

function d2psidr2(r::Float64)
    return (_E0/_b^2)*_td2psidr2(_tr(r))
end

# Effective potential
function psiEff(r::Float64, L::Float64)
    return _E0*_tpsiEff(_tr(r),_tL(L))
end


function dpsiEffdr(r::Float64, L::Float64)
    return (_E0/_b)*_tdpsiEffdr(_tr(r),_tL(L))
end


function d2psiEffdr2(r::Float64, L::Float64)
    return (_E0/_b^2)*_td2psiEffdr2(_tr(r),_tL(L))
end

function _trho(tr::Float64)
    return 1/(1+tr^2)^(5/2)
end

function rho(r::Float64)
    return _rho0*_trho(_tr(r))
end

##################################################
# Compute Lc(E) the angular momentum per unit mass of a circular orbit
##################################################

function _tEta(s::Float64, tE::Float64)
    return (s^2-1)*(1/s - tE)
end

function _sc(tE::Float64)
    t1 = 1/(6*tE)
    t2 = (1+54*tE^2)/(216*tE^3)
    t3 = (1/(4*tE))*sqrt(1+1/(27*tE^2))
    return t1+cbrt(t2+t3)+cbrt(t2-t3)
end

# circular radius
function _rc(tE::Float64)
    sc = _sc(tE)
    return _b*sqrt(sc^2-1)
end

# Circular angular momentum
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

function energy(r::Float64, vr::Float64, vt::Float64)
    return psi(r) + vr^2/(2.0) + vt^2/(2.0)
end

function momentum(r::Float64, vr::Float64, vt::Float64)
    return r*vt
end

##################################################
# Convert (E,L) to (vr,vt) at a given radius r
# Convention vr >= 0
##################################################

function radialVelocitySq(r::Float64,E::Float64, L::Float64)
    return 2*(E-psiEff(r,L))
end

function radialVelocity(r::Float64,E::Float64, L::Float64)
    return sqrt(2*abs((E-psiEff(r,L))))
end

function tangentVelocity(r::Float64,E::Float64, L::Float64)
    return L/r
end

function velocity(r::Float64,E::Float64, L::Float64)
    return sqrt(2*abs((E-psi(r,L))))
end

##################################################
# Solve the degree-3 polynomial equation aX^3 + bX^2 + cX + d = 0
# Gives the real roots
##################################################

function solveRealPoly3(a::Float64, b::Float64, c::Float64, d::Float64)
    p = (3*a*c - b^2)/(3*a^2)
    q = (2*b^3 - 9*a*b*c + 27*a^2*d)/(27*a^3)
    disc = -(4*p^3 + 27*q^2)
    if (disc < 0.0)
        rt1 = cbrt(-q/2 + sqrt((q^2)/4 + (p^3)/27)) + cbrt(-q/2 - sqrt((q^2)/4 + (p^3)/27))
        rt1 -= b/(3*a)
        nbRoot = 1
        return nbRoot, rt1
    else
        rt1 = 2*sqrt(-p/3)*cos((1/3)*acos((3*q)/(2*p) * sqrt(-3/p)))
        rt2 = 2*sqrt(-p/3)*cos((1/3)*acos((3*q)/(2*p) * sqrt(-3/p)) - 2*PI/3)
        rt3 = 2*sqrt(-p/3)*cos((1/3)*acos((3*q)/(2*p) * sqrt(-3/p)) - 4*PI/3)
        rt1 -= b/(3*a)
        rt2 -= b/(3*a)
        rt3 -= b/(3*a)
        nbRoot = 3
        return nbRoot, rt1, rt2, rt3
    end
end

function sign(x::Float64)
    if (x > 0.0)
        return 1.0
    elseif (x < 0.0)
        return -1.0
    else
        return 0.0
    end
end

##################################################
# Bissection Algorithm to find zeroes of functions
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

# Computes numerically rc(L)
# Use bissection
# circular radius
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

# Energy of circular orbit
function Ec(L::Float64)
    rc = _rc_L(L)
    return psiEff(rc,L)
end
