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
