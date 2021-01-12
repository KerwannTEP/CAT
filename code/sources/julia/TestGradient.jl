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
                       (E*(7/2-q)*_H(0,q/2,9/2-q,1,x)
                       - (L^2)/2 * dHdx(0,q/2,9/2-q,1,x))
            end
        end
    end
end

function dfadvrNew(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                   th::Float64, phi::Float64, q::Float64)

    Ea = _Ea(r,vr,vt,vp,th,phi)
    La = _La(r,vr,vt,vp,th,phi)
    return (-vr+vp*cos(th))*dFdE(Ea,La,q)
end

function dfadvrNumRelErrorNew(r::Float64, vr::Float64, vt::Float64, vp::Float64, 
                              th::Float64, phi::Float64, q::Float64, 
                              eps::Float64=0.000001)

    Ea1 = _Ea(r,vr-eps,vt,vp,th,phi)
    La1 = _La(r,vr-eps,vt,vp,th,phi)

    Ea2 = _Ea(r,vr+eps,vt,vp,th,phi)
    La2 = _La(r,vr+eps,vt,vp,th,phi)

    num1 = DistFunction(Ea1,La1,q)
    num2 = DistFunction(Ea2,La2,q)

    real = dfadvrNew(r,vr,vt,vp,th,phi,q)
    return abs(((num2-num1)/(2*eps)-real)/real)
end
