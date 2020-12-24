##################################################
# Computation of the local velocity deflections
##################################################

# Define a Coulomb logarithm
# logCoulomb

function _I1(r::Float64, vr::Float64, vt::Float64)
# integral I1 of the notes
end

function _I2(r::Float64, vr::Float64, vt::Float64)
# integral I2 of the notes
end

function localVelocityChange(r::Float64, vr::Float64, vt::Float64,
                             m_field::Float64, m_test::Float64=0.0)

    cst = 4.0*pi*G^2*m_field
    dvPar  = cst*(m_field+m_test)*logCoulomb
    dvPar2 = cst* m_field 
    dvTan2 = cst* m_field        *(2.0*logCoulomb-1.0)
    I1 = _I1(r,vr,vt)
    I2 = _I2(r,vr,vt)
    dvPar  *= I1
    dvPar2 *= I2
    dvTan2 *= I2
    return dvPar, dvPar2, dvTan2
end