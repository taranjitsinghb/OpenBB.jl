# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:25:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T19:09:39+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# linear constraints set
mutable struct LinearConstraintSet{T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}} <: AbstractConstraintSet
    A::T
    loBs::Array{Float64,1}
    upBs::Array{Float64,1}
    sosIndices::Array{Int,1}
end

function LinearConstraintSet(;A::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},
                    loBs::Array{Float64,1},
                    upBs::Array{Float64,1},
                    sosIndices::Array{Int,1}=Int[])::LinearConstraintSet
    return LinearConstraintSet(A,loBs,upBs,sosIndices)
end


function get_numConstraints(constraintSet::LinearConstraintSet)::Int
    return size(constraintSet.A,1)
end

import SparseArrays.sparse
function sparse(constraintSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(sparse(constraintSet.A),constraintSet.loBs,constraintSet.upBs,constraintSet.sosIndices)
end

function get_sparsity(constraintSet::LinearConstraintSet{SparseMatrixCSC{Float64,Int}})::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(constraintSet.A)[1:2]
end

function get_sparsity(constraintSet::LinearConstraintSet{SparseMatrixCSC{Float64,Int}},index::Int)::Array{Int,1}
    return findnz(constraintSet.A[index,:])[1]
end


import Base.permute!
function permute!(constraintSet::LinearConstraintSet{T},indices::Array{Int,1})::Nothing where T <: Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    constraintSet.A = constraintSet.A[indices,:]
    permute!(constraintSet.loBs,indices)
    permute!(constraintSet.upBs,indices)
    return
end


import Base.append!
function append!(constraintSet1::LinearConstraintSet,constraintSet2::LinearConstraintSet)::Nothing
        constraintSet1.A = vcat(constraintSet1.A,constraintSet2.A)
        append!(constraintSet1.loBs,constraintSet2.loBs)
        append!(constraintSet1.upBs,constraintSet2.upBs)
    return
end
