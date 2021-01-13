using SpecialFunctions
using HypergeometricFunctions


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
                return 3.0*gamma(6.0-q)/(2.0*(2.0*pi)^(5/2)*gamma(q/2.0)) *
                       E^(7/2-q) * _H(0.0,q/2.0,9/2-q,1.0,L^2/(2.0*E))
            end
        end
    end
end

function psi(r::Float64) 
    return 1.0/sqrt(1.0+r^2)
end


function _Ea(r::Float64, vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    v  = sqrt(vr^2  + vt^2)
    return psi(r) - (1/2)*(v^2+vp^2-2*vp*(vr*cos(th)+vt*sin(th)*cos(phi)))
end

function _La(r::Float64, vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    return r*sqrt(vt^2+vp^2*sin(th)^2-2*vt*vp*sin(th)*cos(phi))
end

function dfadvr(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vr+vp*cos(th))*dFdE(Ea,La,q)
end

function dfadvrNumRelError(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64, eps::Float64=0.000001)

    Ea1 = _Ea(r,vr-eps,vt,vp,th,phi)
    La1 = _La(r,vr-eps,vt,vp,th,phi)

    Ea2 = _Ea(r,vr+eps,vt,vp,th,phi)
    La2 = _La(r,vr+eps,vt,vp,th,phi)

    num1 = DistFunction(Ea1,La1,q)
    num2 = DistFunction(Ea2,La2,q)

    real = dfadvr(r,vr,vt,vp,th,phi,q)
    return abs(((num2-num1)/(2*eps)-real)/real)
end

function dfadvt(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vt+vp*sin(th)*cos(phi))*(dFdE(Ea,La,q)-r/La *dFdL(Ea,La,q))
end

function dfadvtNumRelError(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64, eps::Float64=0.000001)

    Ea1 = _Ea(r,vr,vt-eps,vp,th,phi)
    La1 = _La(r,vr,vt-eps,vp,th,phi)

    Ea2 = _Ea(r,vr,vt+eps,vp,th,phi)
    La2 = _La(r,vr,vt+eps,vp,th,phi)

    num1 = DistFunction(Ea1,La1,q)
    num2 = DistFunction(Ea2,La2,q)

    real = dfadvt(r,vr,vt,vp,th,phi,q)
    return abs(((num2-num1)/(2*eps)-real)/real)
end

function d2fadvr2(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                  th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return -dFdE(Ea,La,q) + (-vr+vp*cos(th))^2*d2FdE2(Ea,La,q)
end

function d2fadvr2NumRelError(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64, eps::Float64=0.001)

    Ea1 = _Ea(r,vr-eps,vt,vp,th,phi)
    La1 = _La(r,vr-eps,vt,vp,th,phi)

    Ea2 = _Ea(r,vr+eps,vt,vp,th,phi)
    La2 = _La(r,vr+eps,vt,vp,th,phi)

    Ea0 = _Ea(r,vr,vt,vp,th,phi)
    La0 = _La(r,vr,vt,vp,th,phi)

    num1 = DistFunction(Ea1,La1,q)
    num2 = DistFunction(Ea2,La2,q)
    num0 = DistFunction(Ea0,La0,q)

    real = d2fadvr2(r,vr,vt,vp,th,phi,q)
    return abs(((num2+num1-2*num0)/(eps^2)-real)/real)
end

function d2fadvrdvt(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                    th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return ((-vr+vp*cos(th))*(-vt+vp*sin(th)*cos(phi))
           *(d2FdE2(Ea,La,q) - r/La *d2FdEdL(Ea,La,q)))
end

function d2fadvrdvtNumRelError(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64, eps::Float64=0.001)


    num1 = dfadvr(r,vr,vt-eps,vp,th,phi,q)
    num2 = dfadvr(r,vr,vt+eps,vp,th,phi,q)

    real = d2fadvrdvt(r,vr,vt,vp,th,phi,q)

    return abs(((num2-num1)/(2*eps)-real)/real)
end

function d2fadvt2(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                  th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-dFdE(Ea,La,q) + r/La *dFdL(Ea,La,q) +(-vt+vp*sin(th)*cos(phi))^2
           *(d2FdE2(Ea,La,q) - 2*r/La *d2FdEdL(Ea,La,q)
             -r^2/La^3 *dFdL(Ea,La,q) + r^2/La^2 *d2FdL2(Ea,La,q)))
end

function d2fadvt2NumRelError(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64, eps::Float64=0.001)

    Ea1 = _Ea(r,vr,vt-eps,vp,th,phi)
    La1 = _La(r,vr,vt-eps,vp,th,phi)

    Ea2 = _Ea(r,vr,vt+eps,vp,th,phi)
    La2 = _La(r,vr,vt+eps,vp,th,phi)

    Ea0 = _Ea(r,vr,vt,vp,th,phi)
    La0 = _La(r,vr,vt,vp,th,phi)

    num1 = DistFunction(Ea1,La1,q)
    num2 = DistFunction(Ea2,La2,q)
    num0 = DistFunction(Ea0,La0,q)

    real = d2fadvt2(r,vr,vt,vp,th,phi,q)
    return abs(((num2+num1-2*num0)/(eps^2)-real)/real)
end



function dEadvr(r::Float64, vr::Float64, vt::Float64, 
             vp::Float64, th::Float64, phi::Float64)

    return -vr + vp*cos(th)
end

function dEadvrNum(r::Float64, vr::Float64, vt::Float64, 
                   vp::Float64, th::Float64, phi::Float64,
                   eps::Float64=0.000001)

    num1 = _Ea(r,vr-eps,vt,vp,th,phi)
    num2 = _Ea(r,vr+eps,vt,vp,th,phi)
    return (num2-num1)/(2*eps)
end

function _H(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64)
    if (x <= 1)
        return gamma(a+b)/(gamma(c-a)*gamma(a+d)) * x^a *
               _₂F₁(a+b,1.0+a-c,a+d,x)
    else
        return gamma(a+b)/(gamma(d-b)*gamma(b+c)) * 1.0/x^b *
               _₂F₁(a+b,1.0+b-d,b+c,1.0/x)
    end
end

# problem with q!=0

function dHdx(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64)
    if (x<=1)
        return gamma(a+b)*x^(a-1)/(gamma(c-a)*gamma(a+d)) *
               ( a*_₂F₁(a+b,1+a-c,a+d,x)
               + ((a+b)*(1+a-c))/(a+d) *x*_₂F₁(a+b+1,a-c+2,a+d+1,x))
    else
        return gamma(a+b)*x^(-b-1)/(gamma(d-b)*gamma(b+c)) *
               ( (-b)*_₂F₁(a+b,1+b-d,b+c,1/x)
               - ((a+b)*(1+b-d))/(b+c) *(1/x)*_₂F₁(a+b+1,b-d+2,b+c+1,1/x))
    end
end

function dHdxNum(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64, eps::Float64=0.00001)
    return (_H(a,b,c,d,x+eps)-_H(a,b,c,d,x-eps))/(2*eps)
end

function d2Hdx2(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64)
    if (x<=1)
        return gamma(a+b)*x^(a-2)/(gamma(c-a)*gamma(a+d)) *
               ( a*(a-1)*_₂F₁(a+b,1+a-c,a+d,x)
               + (2*a)*(a+b)*(1+a-c)/(a+d) *x*_₂F₁(a+b+1,a-c+2,a+d+1,x)
               + (a+b)*(1+a-c)/(a+d)*(a+b+1)*(2+a-c)/(a+d+1) *x^2*_₂F₁(a+b+2,a-c+3,a+d+2,x))
    else
        return gamma(a+b)*x^(-b-2)/(gamma(d-b)*gamma(b+c)) *
               ( b*(b+1)*_₂F₁(a+b,1+b-d,b+c,1/x)
               + (2*b+2)*(a+b)*(1+b-d)/(b+c) *(1/x)*_₂F₁(a+b+1,b-d+2,b+c+1,1/x)
               + (a+b)*(1+b-d)/(b+c)*(a+b+1)*(2+b-d)/(b+c+1) *(1/x^2)*_₂F₁(a+b+2,b-d+3,b+c+2,1/x))
    end
end

function d2Hdx2Num(a::Float64, b::Float64, c::Float64, d::Float64, x::Float64, eps::Float64=0.001)
    return (_H(a,b,c,d,x+eps)+_H(a,b,c,d,x-eps)-2*_H(a,b,c,d,x))/(eps^2)
end

function dFdE(E::Float64, L::Float64, q::Float64)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else
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
                return 3*gamma(6-q)*E^(3/2-q)/(2*(2*pi)^(5/2)*gamma(q/2))*
                       (E*(7/2-q)*_H(0.0,q/2,9/2-q,1.0,x)
                       - (L^2)/2 * dHdx(0.0,q/2,9/2-q,1.0,x))
            end
        end
    end
end

function d2FdE2(E::Float64, L::Float64, q::Float64)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else
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
                return 3*gamma(6-q)*E^(3/2-q)/(2*(2*pi)^(5/2)*gamma(q/2))*
                       ((7/2-q)*(5/2-q)*_H(0.0,q/2,9/2-q,1.0,x)
                       - (5-2*q)*(L^2)/(2*E) * dHdx(0.0,q/2,9/2-q,1.0,x)
                       + (L^4)/(4*E^2) * d2Hdx2(0.0,q/2,9/2-q,1.0,x))
            end
        end
    end
end

function d2FdEdL(E::Float64, L::Float64, q::Float64)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else
        x = L^2/(2*E)
        if (q == 0.0) # Isotropic case
            return 0.0
        else
            if (q == 2.0)
                if (x <= 1.0)
                    return -18*L/(2*pi)^3 * (2*E-L^2)^(-1/2)
                else
                    return 0.0
                end
            else
                return 3*gamma(6-q)*L*E^(1/2-q)/(2*(2*pi)^(5/2)*gamma(q/2))*
                       ((5/2-q)*E * dHdx(0.0,q/2,9/2-q,1.0,x)
                        -(L^2)/2* d2Hdx2(0.0,q/2,9/2-q,1.0,x))
            end
        end
    end
end

function dFdL(E::Float64, L::Float64, q::Float64)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else
        x = L^2/(2*E)
        if (q == 0.0) # Isotropic case
            return 0.0
        else
            if (q == 2.0)
                if (x <= 1.0)
                    return -18*L/(2*pi)^3 * (2*E-L^2)^(1/2)
                else
                    return 0.0
                end
            else
                return 3*gamma(6-q)*E^(5/2-q)/(2*(2*pi)^(5/2)*gamma(q/2))*
                       (L * dHdx(0.0,q/2,9/2-q,1.0,x))
            end
        end
    end
end

function d2FdL2(E::Float64, L::Float64, q::Float64)
    if (E <= 0.0 || L <= 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else
        x = L^2/(2*E)
        if (q == 0.0) # Isotropic case
            return 0.0
        else
            if (q == 2.0)
                if (x <= 1.0)
                    return (-18/(2*pi)^3*(2*E-L^2)^(1/2)
                           +18*L^2/(2*pi)^3*(2*E-L^2)^(-1/2))
                else
                    return 0.0
                end
            else
                return 3*gamma(6-q)*E^(3/2-q)/(2*(2*pi)^(5/2)*gamma(q/2))*
                       (E * dHdx(0.0,q/2,9/2-q,1.0,x)
                        +L^2 * d2Hdx2(0.0,q/2,9/2-q,1.0,x))
            end
        end
    end
end

function dFdLNum(E::Float64, L::Float64, q::Float64, eps::Float64=0.001)
    num1 = DistFunction(E,L-eps,q)
    num2 = DistFunction(E,L+eps,q)
    println(num1)
    println(num2)
    return (num2-num1)/(2*eps)
end

function dFdENum(E::Float64, L::Float64, q::Float64, eps::Float64=0.001)
    num1 = DistFunction(E-eps,L,q)
    num2 = DistFunction(E+eps,L,q)
    println(num1)
    println(num2)
    return (num2-num1)/(2*eps)
end

function d2FdE2Num(E::Float64, L::Float64, q::Float64, eps::Float64=0.001)
    num1 = DistFunction(E-eps,L,q)
    num2 = DistFunction(E+eps,L,q)
    num0 = DistFunction(E,L,q)
    return (num1+num2-2*num0)/(eps^2)
end

function d2FdEdLNum(E::Float64, L::Float64, q::Float64, eps::Float64=0.001)
    num1 = dFdL(E-eps,L,q)
    num2 = dFdL(E+eps,L,q)
    return (num2-num1)/(2*eps)
end

function d2FdL2Num(E::Float64, L::Float64, q::Float64, eps::Float64=0.0001)
    num1 = DistFunction(E,L-eps,q)
    num2 = DistFunction(E,L+eps,q)
    num0 = DistFunction(E,L,q)
    return (num1+num2-2*num0)/(eps^2)
end

function d2FdL2Num1(E::Float64, L::Float64, q::Float64, eps::Float64=0.001)
    num1 = dFdL(E,L-eps,q)
    num2 = dFdL(E,L+eps,q)
    return (num2-num1)/(2*eps)
end

function dfadvrNew(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vr+vp*cos(th))*dFdE(Ea,La,q)
end

function dfadvrNumNew(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                      th::Float64, phi::Float64, q::Float64, 
                      eps::Float64=0.001)

    Ea1 = _Ea(r,vr-eps,vt,vp,th,phi)
    La1 = _La(r,vr-eps,vt,vp,th,phi)

    Ea2 = _Ea(r,vr+eps,vt,vp,th,phi)
    La2 = _La(r,vr+eps,vt,vp,th,phi)

    num1 = DistFunction(Ea1,La1,q)
    num2 = DistFunction(Ea2,La2,q)

#    real = dfadvrNew(r,vr,vt,vp,th,phi,q)
    return (num2-num1)/(2*eps)
end

#problem dans les definition?
