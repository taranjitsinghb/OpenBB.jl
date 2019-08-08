# structure used to hold the settings for OSQP
mutable struct MOSEKsettings <: AbstractSettings
LOG::Int
end


function MOSEKsettings(; LOG=1)::MOSEKsettings

    return MOSEKsettings(LOG)
end

mutable struct MOSEKworkspace <: AbstractWorkspace
    #objectives
    Q::SparseMatrixCSC{Float64}
    L::Array{Float64,1}

    #constraints
    A::SparseMatrixCSC{Float64}
    cnsLoBs::Array{Float64,1}
    cnsUpBs::Array{Float64,1}

    #variables
    varLoBs::Array{Float64,1}
    varUpBs::Array{Float64,1}
    #X::Array{Float64,1}

    #workspace
    model::JUMP_MOSEK.Model
    settings::MOSEKsettings
end
