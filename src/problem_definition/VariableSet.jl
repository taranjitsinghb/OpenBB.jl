# @Author: Massimo De Mauri <massimo>
# @Date:   2019-04-23T15:05:12+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: VariableSet.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-17T15:47:54+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractVariableSet end
struct NullVariableSet <: AbstractVariableSet end



# struct to store the variables data
struct VariableSet <: AbstractVariableSet
    loBs::Array{Float64,1}
    upBs::Array{Float64,1}
    val::Array{Float64,1}
    dscIndices::Array{Int,1}
    sos1Groups::Array{Int,1} # assume group -1 as no group
    pseudoCosts::Array{Float64,2}
end

function VariableSet(;loBs::Array{Float64,1},upBs::Array{Float64,1},val::Array{Float64,1}=Float64[],
                      dscIndices::Array{Int,1}=Int[],sos1Groups::Array{Int,1}=Int[],pseudoCosts=Float64[])::VariableSet

    # check the correctness of inputs
    if length(sos1Groups) == 0
        sos1Groups = repeat([-1],length(dscIndices))
    elseif length(sos1Groups) != length(dscIndices)
        @error "sos1Groups should either be empty or have the same length of dscIndices"
    end

    if size(pseudoCosts,1) == 0
        pseudoCosts = 1e-4*ones(length(dscIndices),2)
    elseif size(pseudoCosts,1) != length(dscIndices) || size(pseudoCosts,2) != 2
        @error "pseudoCosts should either be empty or have size: length(dscIndices) x 2 "
    end

    if length(val) == 0
        @assert length(loBs) == length(upBs)
        val = Array{Float64,1}(undef,length(loBs))
        for i in 1:length(loBs)
            scenario = (loBs[i]>-Inf) + 2*(upBs[i]<Inf)
            if scenario == 3
                val[i] = .5*(loBs[i] + upBs[i])
            elseif scenario == 2
                val[i] = upBs[i]
            elseif scenario == 1
                val[i] = loBs[i]
            else
                val[i] = 0
            end
        end
    else
        @assert length(loBs) == length(upBs) == length(val)
    end

    return VariableSet(loBs,upBs,val,dscIndices,sos1Groups,pseudoCosts)
end

function EmptyVarSet()::VariableSet
    return VariableSet(Float64[],Float64[],Float64[],Int[],Int[],Array{Float64,2}(undef,0,0))
end



function get_numVariables(variableSet::VariableSet)::Int
    return length(variableSet.loBs)
end

function get_numDiscreteVariables(variableSet::VariableSet)::Int
    return length(variableSet.dscIndices)
end
