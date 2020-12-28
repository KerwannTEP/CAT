using HypergeometricFunctions
using SpecialFunctions

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