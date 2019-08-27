# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:33:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T14:29:10+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# constructors and copy functions (Fundamental. These are used in Branch and Bound)
# named constructor
function LinearObjective(;L::T)::LinearObjective where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return LinearObjective(L)
end

# type conversions
function LinearObjective{T}(objective::LinearObjective{T})::LinearObjective{T} where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return objective
end

function LinearObjective{T1}(objective::QuadraticObjective{T1,T2})::LinearObjective{T1} where T1<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}} where T2<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert all(objective.Q .== 0)
    return LinearObjective(objective.L)
end

import Base.copy
function copy(objective::LinearObjective)::LinearObjective
    return LinearObjective(objective.L)
end

import Base.deepcopy
function deepcopy(objective::LinearObjective)::LinearObjective
    return LinearObjective(deepcopy(objective.L))
end

import SparseArrays.sparse
function sparse(objective::LinearObjective)::LinearObjective
    return LinearObjective(sparse(objective.L))
end

# inspect functions (Fundamental. These are used in Branch and Bound)
function get_numVariables(objective::LinearObjective)::Int
    return size(objective.L,1)
end

function get_sparsity(objective::LinearObjective)::Array{Int,1}
    return findnz(objective.L)[1]
end


# update functions (Not Fundamental. These are used only during problem update)
function insert_variables!(objective::LinearObjective,numVariables::Int,insertionPoint::Int)::Nothing
    splice!(objective.L,insertionPoint:insertionPoint-1,zeros(numVariables,1))
    return
end

function append_variables!(objective::LinearObjective,numVariables::Int)::Nothing
    insert_variables!(objective,numVariables,get_numVariables(objective)+1)
    return
end


function remove_variables!(objective::LinearObjective,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objective)))
    objective.L = objective.L[toKeep]
    return
end


import Base.+
function +(objective1::LinearObjective,objective2::T)::LinearObjective where T<:AbstractObjective

    if objective2 isa NullObjective
        return objective1
    else
        @assert get_numVariables(objective1) == get_numVariables(objective2)
        objective2 = LinearObjective(objective2)
        return LinearObjective(objective1.L+objective2.L)
    end
end
function +(objective1::T,objective2::LinearObjective)::LinearObjective where T<:AbstractObjective

    if objective1 isa NullObjective
        return objective2
    else
        @assert get_numVariables(objective1) == get_numVariables(objective2)
        objective1 = LinearObjective(objective1)
        return LinearObjective(objective1.L+objective2.L)
    end
end
