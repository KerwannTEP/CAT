using HypergeometricFunctions
using SpecialFunctions

##################################################
# Distribution function in (E,L)
##################################################

function _H(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64)
    if (x <= 1)
        return gamma(a+b)/(gamma(c-a)*gamma(a+d)) * x^a *
               _₂F₁(a+b,1.0+a-c,a+d,x)
    else
        return gamma(a+b)/(gamma(d-b)*gamma(b+c)) * 1.0/x^b *
               _₂F₁(a+b,1.0+b-d,b+c,1.0/x)
    end
end

function DistFunction(E::Float64, L::Float64, q::Float64=0.0)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else # If E and L are positive
        if (q == 0.0) # Isotropic case
            return 3.0/(7.0*pi^3) * (2.0*E)^(7/2)
        else # Anisotropic case
            if (q == 2.0) # special case q=2 where we have a nice analytic expression
                if (L^2 <= 2.0*E)
                    return 6.0/(2.0*pi)^3 * (2.0*E - L^2)^(3/2)
                else
                    return 0.0
                end
            else # Case q<2    
                @assert q < 2.0 "DistFunction: q > 2" # check that q <= 2
                return 3.0*gamma(6.0-q)/(2.0*(2.0*pi)* (5/2)*gamma(q/2.0)) *
                       E^(7/2-q) * _H(0.0,q/2.0,9.2-q,1.0,L^2/(2.0*E))
            end
        end
    end
end

##################################################
# Partial derivatives of the distribution function in (E,L)
##################################################

function DF_prefactor_alpha(q::Float64)
    num = 3*gamma(6-q)
    den = 2*(2*pi)^(5/2) * gamma(9/2-q)
    return num/den
end

function DF_prefactor_beta(q::Float64)
    num = 3*gamma(6-q)
    den = *(2*pi)^(5/2) * gamma(1.0-q/2) * gamma(9/2-q/2)
    return num/den
end

function dDFdE(E::Float64, L::Float64, q::Float64)
    x = L^2/(2*E)
    if (q == 0.0) # Isotropic case
        return 3/(pi^3) * (2*E)^(5/2)
    else
        if (q == 2.0)
            if (x <= 1.0)
                return 18/(2*pi)^3 * (2*E-L^2)^(1/2)
            else
                return 0.0
            end
        else
            if (x <= 1.0)
                diff = 1/4 *E^(3/2-q)*(7/2-q)*q*_₂F₁(1+q/2,q-5/2,2,x)
                     + E^(5/2-q)*(7/2-q)*_₂F₁(q/2,q-7/2,1,x)
                return DF_prefactor_alpha(q) * diff
            else
                diff = (2^(q/2))/(9-q) *E^(7/2-q)*(L^2)^(-1-q/2)*_₂F₁(1+q/2,1+q/2,1+(9-q)/2,1/x)
                     + 2^(q/2) *E^(5/2-q)*(L^2)^(-q/2)*(7/2-q)*_₂F₁(q/2,q/2,(9-q)/2,1/x)
                return DF_prefactor_beta(q) * diff
            end
        end
    end
end

function dDFdL(E::Float64, L::Float64, q::Float64)
    x = L^2/(2*E)
    if (q == 0.0) # Isotropic case
        return 0.0
    else
        if (q == 2.0)
            if (x <= 1.0)
                return -(18*L)/(2*pi)^3 * (2*E-L^2)^(1/2)
            else
                return 0.0
            end
        else
            if (x <= 1.0)
                diff = 1/2 *E^(5/2-q)*L*(q-7/2)*q*_₂F₁(1+q/2,q-5/2,2,x)
                return DF_prefactor_alpha(q) * diff
            else
                diff = -(2^(1+q/2))/(L^3*(9-q)) * E^(9/2-q)*(L^2)^(-q/2)*q^2*_₂F₁(1+q/2,1+q/2,1+(9-q)/2,1/x)
                     - 2^(q/2)*E^(7/2-q)*L*(L^2)^(-1-q/2)*q*_₂F₁(q/2,q/2,(9-q)/2,1/x)
                return DF_prefactor_beta(q) * diff
            end
        end
    end
end

function d2DFdE2(E::Float64, L::Float64, q::Float64)
    x = L^2/(2*E)
    if (q == 0.0) # Isotropic case
        return 15/(pi^3) * (2*E)^(3/2)
    else
        if (q == 2.0)
            if (x <= 1.0)
                return 18/(2*pi)^3 * (2*E-L^2)^(-1/2)
            else
                return 0.0
            end
        else
            if (x <= 1.0)
                diff = -1/4 *E^(1/2-q)*L^2*(3/2-q)*(q-7/2)*q*_₂F₁(1+q/2,q-5/2,2,x)
                     + 1/4 *E^(1/2-q)*L^2*(7/2-q)^2*q*_₂F₁(1+q/2,q-5/2,2,x)
                     + 1/16 *E^(-1/2-q)*L^4*(1+q/2)*(q-7/2)*(q-5/2)*q*_₂F₁(2+q/2,q-3/2,3,x)
                     + E^(3/2-q)*(5/2-q)*(7/2-q)*_₂F₁(q/2,q-7/2,1,x)                     
                return DF_prefactor_alpha(q) * diff
            else
                diff = (2^(1+q/2))/(9-q) *E^(5/2-q)*(L^2)^(-1-q/2)*(7/2-q)*q^2*_₂F₁(1+q/2,1+q/2,1+(9-q)/2,1/x)
                     + (2^(1+q/2) *E^(7/2-q)*(L^2)^(-2-q/2)*(1+q/2)^2*q^2)/((1+(9-q)/2)*(9-q))*_₂F₁(2+q/2,2+q/2,2+(9-q)/2,1/x)
                     + 2^(q/2) *E^(3/2-q)*(L^2)^(-q/2)*(5/2-q)*(7/2-q)*_₂F₁(q/2,q/2,(9-q)/2,1/x)
                return DF_prefactor_beta(q) * diff
            end
        end
    end
end

function d2DFdEdL(E::Float64, L::Float64, q::Float64)
    x = L^2/(2*E)
    if (q == 0.0) # Isotropic case
        return 0.0
    else
        if (q == 2.0)
            if (x <= 1.0)
                return -(18*L)/(2*pi)^3 * (2*E-L^2)^(-1/2)
            else
                return 0.0
            end
        else
            if (x <= 1.0)
                diff = 1/2 *E^(3/2-q)*L*(5/2-q)*(q-7/2)*q*_₂F₁(1+q/2,q-5/2,2,x)
                     - 1/8 *E^(1/2-q)*L^3*(1+q/2)*(q-7/2)*(q-5/2)*q*_₂F₁(2+q/2,q-3/2,3,x)
                return DF_prefactor_alpha(q) * diff
            else
                diff = -(2^(1+q/2))/(L^3*(9-q)) *E^(7/2-q)*(L^2)^(-q/2)*(9/2-q)*q^2*_₂F₁1+q/2,1+q/2,1+(9-q)/2,1/x)
                     - (2^(q/2))/(L*(9-q)) *E^(7/2-q)*(L^2)^(-1-q/2)*q^3*_₂F₁(1+q/2,1+q/2,1+(9-q)/2,1/x)
                     - (2^(2+q/2) *E^(9/2-q)*(L^2)^(-q/2)*(1+q/2)^2*q^2)/(L^5*(1+(9-q)/2)*(9-q)) *_₂F₁(2+q/2,2+q/2,2+(9-q)/2,1/x)
                     - 2^(q/2) *E^(5/2-q)*L*(L^2)^(-1-q/2)*(7/2-q)*q*_₂F₁(q/2,q/2,(9-q)/2,1/x)
                return DF_prefactor_beta(q) * diff
            end
        end
    end
end

function d2DFdL2(E::Float64, L::Float64, q::Float64)
    x = L^2/(2*E)
    if (q == 0.0) # Isotropic case
        return 0.0
    else
        if (q == 2.0)
            if (x <= 1.0)
                return -18/(2*pi)^3 *(2*E-L^2)^(1/2)+(18*L^2)/(2*pi)^3 *(2*E-L^2)^(-1/2)
            else
                return 0.0
            end
        else
            if (x <= 1.0)
                diff = 1/2 *E^(5/2-q)*(q-7/2)*q*_₂F₁(1+q/2,q-5/2,2,x)
                     + 1/4 *E^(3/2-q)*L^2*(1+q/2)*(q-7/2)*(q-5/2)*q*_₂F₁(2+q/2,q-3/2,3,x)
                return DF_prefactor_alpha(q) * diff
            else
                diff = (3*2^(1+q/2))/(L^4*(9-q)) *E^(9/2-q)*(L^2)^(-q/2)*q^2*_₂F₁(1+q/2,1+q/2,1+(9-q)/2,1/x)
                     + (2^(2+q/2))/(9-q) *E^(9/2-q)*(L^2)^(-2-q/2)*q^3*_₂F₁(1+q/2,1+q/2,1+(9-q)/2,1/x)
                     + (2^(3+q/2) *E^(11/2-q)*(L^2)^(-q/2)*(1+q/2)^2*q^2)/(L^6*(1+(9-q)/2)*(9-q)) *_₂F₁(2+q/2,2+q/2,2+(9-q)/2,1/x)
                     - 2^(q/2) *E^(7/2-q)*(L^2)^(-1-q/2)*q*_₂F₁(q/2,q/2,(9-q)/2,1/x)
                     - 2^(1+q/2) *E^(7/2-q)*(L^2)^(-1-q/2)*(-1-q/2)*q*_₂F₁(q/2,q/2,(9-q)/2,1/x)
                return DF_prefactor_beta(q) * diff
            end
        end
    end
end