# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:33:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-25T17:09:38+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# constructors and copy functions (Fundamental. These are used in Branch and Bound)
# named constructor
function LinearObjective(;L::T)::LinearObjective{T} where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return LinearObjective(L)
end

function LinearObjective(;L::T)::LinearObjective{T} where T<:Union{Array{Float64,2},SparseVector{Float64,Int}}
    return LinearObjective(L)
end

# type conversions
function LinearObjective{T}(objective::LinearObjective{T})::LinearObjective{T} where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return objective
end

function LinearObjective{T}(objective::LinearObjective{T})::LinearObjective{T} where T<:Union{Array{Float64,2},SparseVector{Float64,Int}}
    return objective
end

function LinearObjective{T}(objective::LinearObjective)::LinearObjective{T} where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    return LinearObjective(T(objective.L))
end


function LinearObjective{T}(objective::QuadraticObjective)::LinearObjective{T} where T<:Union{Array{Float64,1},SparseVector{Float64,Int}}
    @assert all(objective.Q .== 0)
    return LinearObjective(T(objective.L))
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
    return size(objective.L,1)*size(objective.L,2)
end

function get_sparsity(objective::LinearObjective)::Array{Int,1}
    return findnz(sparse(objective.L))[1:2]
end


# update functions (Not Fundamental. These are used only during problem update)
function insert_variables!(objective::LinearObjective,numVariables::Int,insertionPoint::Int)::Nothing
    @assert numVariables>=0
    @assert 0<=insertionPoint<=get_numVariables(objective)+1
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
function +(objective1::LinearObjective,objective2::LinearObjective)::LinearObjective
    @assert get_numVariables(objective1) == get_numVariables(objective2)
    return LinearObjective(objective1.L+objective2.L)
end

function +(objective1::LinearObjective,objective2::T)::LinearObjective where T<:AbstractObjective

    if objective2 isa NullObjective
        return objective1
    else
        return objective1 + LinearObjective(objective2)
    end
end


function add!(objective1::LinearObjective,objective2::LinearObjective)::Nothing
    @assert get_numVariables(objective1) == get_numVariables(objective2)
    @. objective1.L += objective2.L
    return
end

function add!(objective1::LinearObjective,objective2::T)::Nothing where T<:AbstractObjective
    if objective2 isa NullObjective
        return
    else
        add!(objective1,LinearObjective(objective2))
        return
    end
end
