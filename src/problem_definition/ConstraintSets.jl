# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ConstraintSets.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-17T23:22:50+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractConstraintSet end; struct NullConstraintSet <: AbstractConstraintSet end

# linear constraints set
mutable struct LinearCns{T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}} <: AbstractConstraintSet
    A::T
    loBs::Array{Float64,1}
    upBs::Array{Float64,1}
    sosIndices::Array{Int,1}
end

function LinearCns(;A::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},
                    loBs::Array{Float64,1},
                    upBs::Array{Float64,1},
                    sosIndices::Array{Int,1}=Int[])::LinearCns
    return LinearCns(A,loBs,upBs,sosIndices)
end


import Base.size
function get_numConstraints(constraintSet::LinearCns)::Int
    return size(constraintSet.A,1)
end

function get_numVariables(constraintSet::LinearCns)::Int
    return size(constraintSet.A,2)
end

import SparseArrays.sparse
function sparse(constraintSet::LinearCns)::LinearCns
    return LinearCns(sparse(constraintSet.A),constraintSet.loBs,constraintSet.upBs,constraintSet.sosIndices)
end


function get_sparsity(constraintSet::LinearCns{SparseMatrixCSC{Float64,Int}})::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(constraintSet.A)[1:2]
end

function get_sparsity(constraintSet::LinearCns{SparseMatrixCSC{Float64,Int}},index::Int)::Array{Int,1}
    return findnz(constraintSet.A[index,:])[1]
end



import Base.permute!
function permute!(constraintSet::LinearCns{T},indices::Array{Int,1})::Nothing where T <: Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    constraintSet.A = constraintSet.A[indices,:]
    permute!(constraintSet.loBs,indices)
    permute!(constraintSet.upBs,indices)
    return
end


import Base.append!
function append!(constraintSet1::LinearCns,constraintSet2::LinearCns)::Nothing
        constraintSet1.A = vcat(constraintSet1.A,constraintSet2.A)
        append!(constraintSet1.loBs,constraintSet2.loBs)
        append!(constraintSet1.upBs,constraintSet2.upBs)
    return
end
