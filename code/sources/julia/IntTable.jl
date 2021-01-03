mutable struct IntTable
    dhdvr::Base.RefValue{Float64}
    dhdvt::Base.RefValue{Float64}
    dgdvt::Base.RefValue{Float64}	
    d2gdvr2::Base.RefValue{Float64}
    d2gdvrdvt::Base.RefValue{Float64}
    d2gdvt2::Base.RefValue{Float64}
    dvPar::Base.RefValue{Float64}
    dvPar2::Base.RefValue{Float64}
    dvTan2::Base.RefValue{Float64}
    dE::Base.RefValue{Float64}
    dE2::Base.RefValue{Float64}
    dL::Base.RefValue{Float64}
    dL2::Base.RefValue{Float64}
    dEdL::Base.RefValue{Float64}
end

function IntTable_create!()
    return IntTable(Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),
                    Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0),Ref(0.0))
end

function IntTable_init!(Table::IntTable)
    Table.dhdvr = 0.0
    Table.dhdvt = 0.0
    Table.dgdvt	= 0.0
    Table.d2gdvr2 = 0.0
    Table.d2gdvrdvt = 0.0
    Table.d2gdvt2 = 0.0
    Table.dvPar = 0.0
    Table.dvPar2 = 0.0
    Table.dvTan2 = 0.0
    Table.dE = 0.0
    Table.dE2 = 0.0
    Table.dL = 0.0
    Table.dL2 = 0.0
    Table.dEdL = 0.0
end