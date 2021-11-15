using HypergeometricFunctions
using SpecialFunctions
using StaticArrays # To have access to static arrays
using Interpolations # To have access to interpolation functions

##################################################
# Distribution function in (E,L)
##################################################

function _H(x::Float64, q::Float64=qCalc)
    if (x <= 1)
        pref = 1/(GAMMA_ca)
        HG   = H_1(x) #tabHyperGeoInt[1](x)#_₂F₁(q/2,q-3.5,1.0,x)
        return pref*HG
    else
        pref = 1/(GAMMA_db*GAMMA_bc)
        HG   = H_3(1/x) #tabHyperGeoInt[3](1/x)#_₂F₁(q/2,q/2,4.5-q/2,1/x)
        return pref*x^(-q/2)*HG
    end
end

function _tF(tE::Float64, tL::Float64)

    if (tE < 0.0 || tL < 0.0) # If E or L are negative, the DF vanishes
        return 0.0
    else # If E and L are positive

        # FIND A WAY TO ONLY SELECT ONE PART OF THE CODE
        # AT THE COMPILATION
        #

        if (qCalc == 0.0) # Isotropic case
             return 3.0/(7.0*PI^3) * (2.0*tE)^(7/2)
        elseif (qCalc == 2.0) # Anisotropic case
            x = tL^2/(2*tE)
            if (x <= 1)
                return 6.0/(2.0*PI)^3 * (2.0*tE - tL^2)^(3/2)
            else
                return 0.0
            end
        else # Case q<2
            x = tL^2/(2*tE)
            return (3.0*GAMMA_6q/(2.0*(2.0*PI)^(5/2)) *
                       tE^(7/2-qCalc) * _H(x,qCalc))
        end

    end
end

function _F(E::Float64, L::Float64)
    tE = _tE(E)
    tL = _tL(L)
    DF  = _tF(tE,tL)
    return _M*_F0*DF
end


##################################################
# Distribution function in (E,L)
# For flux computation
##################################################

function _HdH(x::Float64, q::Float64=qCalc)
    if (x <= 1)
        pref = 1/(GAMMA_ca*GAMMA_ad)
        HG   = tabHyperGeoInt[1](x) #_₂F₁(q/2,q-3.5,1.0,x)
        HGp  = tabHyperGeoInt[2](x) #_₂F₁(q/2+1,q-2.5,2.0,x)
        H = pref*HG
        dH = pref*x^(-1)*((q/2)*(q-3.5)*x*HGp)
        return H, dH
    else
        pref = 1/(GAMMA_db*GAMMA_bc)
        HG   = tabHyperGeoInt[3](1/x) #_₂F₁(q/2,q/2,4.5-q/2,1/x)
        HGp  = tabHyperGeoInt[4](1/x) #_₂F₁(q/2+1,q/2+1,5.5-q/2,1/x)
        H = pref*x^(-q/2)*HG
        dH = pref*x^(-q/2-1)*((-q/2)*HG - (1/x)*HGp*(q/2)*(q/2)/(4.5-q/2))
        return H, dH
    end
end

# Returns tF, dtFdtE, dtFdtL
function _tFdF(tE::Float64, tL::Float64, q::Float64=qCalc)
    if (tE <= 0.0 || tL <= 0.0) # If E or L are negative, the DF vanishes
        DF = 0.0
        dE = 0.0
        dL = 0.0
        return DF, dE, dL
    else
        x = tL^2/(2*tE)
        if (q == 0.0) # Isotropic case
            pref = 3/(PI^3) * (2*tE)^(5/2)
            DF = pref/7 * (2*tE)
            dE = pref
            dL = 0.0
            return DF, dE, dL
        else
            if (q == 2.0)
                if (x <= 1.0)
                    pref = 18/(2*PI)^3*(2*tE-tL^2)^(1/2)
                    DF = (1/3)*pref* (2*tE - tL^2)
                    dE = pref
                    dL = -pref*tL
                    return DF, dE, dL
                else
                    DF = 0.0
                    dE = 0.0
                    dL = 0.0
                    return DF, dE, dL
                end
            else
                pref = 3*GAMMA_6q/(2*(2*PI)^(5/2))*tE^(5/2-q)
                pref0 = pref*tE
                pref1 = pref
                H, dH = _HdH(x,q)

                DF = pref0 * H
                dE = pref1*((7/2-q)*H - x*dH)
                dL = pref1*(tL*dH)
                return DF, dE, dL
            end
        end
    end
end

function _FdF(E::Float64, L::Float64, q::Float64=qCalc)
    tE = _tE(E)
    tL = _tL(L)
    tDF, tdE, tdL = _tFdF(tE,tL,q)
    DF = _M*_F0*tDF
    dE = (_M*_F0/_E0)*tdE
    dL = (_M*_F0/_L0)*tdL
    return DF, dE, dL
end

##################################################
# Interpolate the hypergeometric function
# HG1   = _₂F₁(qCalc/2,qCalc-3.5,1.0,x)
# HGp1  = _₂F₁(qCalc/2+1,qCalc-2.5,2.0,x)
# HG2   = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,1/x)
# HGp2  = _₂F₁(qCalc/2+1,qCalc/2+1,5.5-qCalc/2,1/x)

# valeur en x=1: https://homepage.tudelft.nl/11r49/documents/wi4006/hyper.pdf equation 4
# 2F1(a,b,c,1) = gamma(c)gamma(c-a-b)/(gamma(c-a)gamma(c-b))

# Valeur en x=1
# HG1   = _₂F₁(qCalc/2,qCalc-3.5,1.0,1) = gamma(4.5-1.5*qCalc)/(gamma(1-qCalc/2)gamma(4.5-qCalc))
# HGp1  = _₂F₁(qCalc/2+1,qCalc-2.5,2.0,1) = gamma(3.5-1.5*qCalc)/(gamma(1-qCalc/2)gamma(4.5-qCalc))
# HG2   = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,1) = gamma(4.5-qCalc/2)gamma(4.5-1.5*qCalc)/(gamma(4.5-qCalc)gamma(4.5-qCalc))
# HGp2  = _₂F₁(qCalc/2+1,qCalc/2+1,5.5-qCalc/2,1) = gamma(5.5-qCalc/2)gamma(3.5-1.5*qCalc)/(gamma(4.5-qCalc)gamma(4.5-qCalc))
##################################################

function getHGInt(nbxInt::Int64=1000)
    xminInt = 0.0
    xmaxInt = 1.0
    #####
    rangexInt = range(xminInt,length=nbxInt,xmaxInt)
    tabxInt = collect(rangexInt)
    tabHG1Int = zeros(Float64,nbxInt)
    tabHGp1Int = zeros(Float64,nbxInt)
    tabHG2Int = zeros(Float64,nbxInt)
    tabHGp2Int = zeros(Float64,nbxInt)

    #####
    for indx=2:nbxInt-1
        xloc = tabxInt[indx]
        hg1loc = _₂F₁(qCalc/2,qCalc-3.5,1.0,xloc)
        hgp1loc = _₂F₁(qCalc/2+1,qCalc-2.5,2.0,xloc)
        hg2loc = _₂F₁(qCalc/2,qCalc/2,4.5-qCalc/2,xloc)
        hgp2loc = _₂F₁(qCalc/2+1,qCalc/2+1,5.5-qCalc/2,xloc)
        tabHG1Int[indx] = hg1loc
        tabHGp1Int[indx] = hgp1loc
        tabHG2Int[indx] = hg2loc
        tabHGp2Int[indx] = hgp2loc
    end
    # x=0
    tabHG1Int[1] = 1.0
    tabHGp1Int[1] = 1.0
    tabHG2Int[1] = 1.0
    tabHGp2Int[1] = 1.0

    # x=1
    tabHG1Int[nbxInt] = gamma(4.5-1.5*qCalc)/(gamma(1-qCalc/2)*gamma(4.5-qCalc))
    tabHGp1Int[nbxInt] = gamma(3.5-1.5*qCalc)/(gamma(1-qCalc/2)*gamma(4.5-qCalc))
    tabHG2Int[nbxInt] = gamma(4.5-qCalc/2)*gamma(4.5-1.5*qCalc)/(gamma(4.5-qCalc)*gamma(4.5-qCalc))
    tabHGp2Int[nbxInt] = gamma(5.5-qCalc/2)*gamma(3.5-1.5*qCalc)/(gamma(4.5-qCalc)*gamma(4.5-qCalc))
    #####
    intHG1 = Interpolations.scale(interpolate(tabHG1Int, BSpline(Cubic(Line(OnGrid())))),rangexInt)
    intHGp1 = Interpolations.scale(interpolate(tabHGp1Int, BSpline(Cubic(Line(OnGrid())))),rangexInt)
    intHG2 = Interpolations.scale(interpolate(tabHG2Int, BSpline(Cubic(Line(OnGrid())))),rangexInt)
    intHGp2 = Interpolations.scale(interpolate(tabHGp2Int, BSpline(Cubic(Line(OnGrid())))),rangexInt)
    #####
    return [intHG1, intHGp1, intHG2, intHGp2]
end


const tabHyperGeoInt = SVector{4}(getHGInt())

const H_1 = tabHyperGeoInt[1]
const H_3 = tabHyperGeoInt[3]
