# structure used to hold the settings for OSQP
mutable struct MOSEKsettings <: AbstractSettings
LOG::Int
end


function MOSEKsettings(; LOG=1)::MOSEKsettings

    return MOSEKsettings(LOG)
end

mutable struct MOSEKworkspace <: AbstractWorkspace
    #objectives
    C::Array{Float64,1}

    #constraints
    A::SparseMatrixCSC{Float64}
    b::Array{Float64,1}

    #variables
    X::Array{Float64,1}

    #workspace
    model::JUMP_MOSEK.Model
    settings::MOSEKsettings
end
