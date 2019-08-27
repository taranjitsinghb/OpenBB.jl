# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:25:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T13:44:17+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# constructirs and copy functions (Fundamental. Those are used in Branch and Bound)
# named constructor
function LinearConstraintSet(;A::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},
                              loBs::Array{Float64,1},
                              upBs::Array{Float64,1})::LinearConstraintSet
    return LinearConstraintSet(A,loBs,upBs)
end


import Base.copy
function copy(constraintSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(constraintSet.A,constraintSet.loBs,constraintSet.upBs)
end

import Base.deepcopy
function deepcopy(constraintSet::LinearConstraintSet)::LinearConstraintSet
    return LinearConstraintSet(deepcopy(constraintSet.A),copy(constraintSet.loBs),copy(constraintSet.upBs))
end

# type conversion
function LinearConstraintSet(constraintSet::LinearConstraintSet)::LinearConstraintSet
    return constraintSet
end

import SparseArrays.sparse
function sparse(constraintSet::LinearConstraintSet)::LinearConstraintSet{SparseMatrixCSC{Float64,Int}}
    return LinearConstraintSet(sparse(constraintSet.A),constraintSet.loBs,constraintSet.upBs)
end


# inspect functions (Fundamental. Those are used in Branch and Bound)
function get_numVariables(constraintSet::LinearConstraintSet)::Int
    return size(constraintSet.A,2)
end

function get_size(constraintSet::LinearConstraintSet)::Int
    return size(constraintSet.A,1)
end

function get_bounds(constraintSet::LinearConstraintSet)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (constraintSet.loBs,constraintSet.upBs)
end


# inspect functions (Not fundamental. Those are used only in updating the problem))
function get_sparsity(constraintSet::LinearConstraintSet)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(constraintSet.A)[1:2]
end

function get_sparsity(constraintSet::LinearConstraintSet,index::Int)::Array{Int,1}
    return findnz(workspace.A[index,:])[1]
end

# update functions (Not fundamental. Those are used only in updating the problem)
function update_bounds!(constraintSet::LinearConstraintSet,loBs::Array{Float64,1},upBs::Array{Float64,1})::Nothing
    @assert length(loBs) == length(upBs) == length(constraintSet.loBs) == length(constraintSet.upBs)
    constraintSet.loBs = loBs
    constraintSet.upBs = upBs
    return
end

function remove_variables!(constraintSet::LinearConstraintSet,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(constraintSet)))
    constraintSet.A = constraintSet.A[:,toKeep]
    return
end

function insert_variables!(constraintSet::LinearConstraintSet,numVariables::Int,insertionPoint::Int)::Nothing
    constraintSet.A = hcat(constraintSet.A[:,1:insertionPoint-1],zeros(size(constraintSet.A,1),numVariables),constraintSet.A[:,insertionPoint:end])
    return
end

function append_variables!(constraintSet::LinearConstraintSet,numVariables::Int)::Nothing
    insert_variables!(constraintSet,numVariables,get_size(constraintSet)+1)
    return
end

function remove_constraints!(constraintSet::LinearConstraintSet,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_size(constraintSet)))
    constraintSet.A = constraintSet.A[toKeep,:]
    constraintSet.loBs = constraintSet.loBs[toKeep]
    constraintSet.upBs = constraintSet.upBs[toKeep]
    return
end

function insert_constraints!(constraintSet1::LinearConstraintSet,constraintSet2::LinearConstraintSet,insertionPoint::Int)::Nothing
        constraintSet1.A = vcat(constraintSet1.A[1:insertionPoint-1,:],constraintSet2.A,constraintSet1.A[insertionPoint:end,:])
        splice!(constraintSet1.loBs,insertionPoint:insertionPoint-1,copy(constraintSet2.loBs))
        splice!(constraintSet1.upBs,insertionPoint:insertionPoint-1,copy(constraintSet2.upBs))
    return
end

function append_constraints!(constraintSet1::LinearConstraintSet,constraintSet2::LinearConstraintSet)::Nothing
    insert_constraints!(constraintSet1,constraintSet2,get_size(constraintSet1)+1)
    return
end

function permute_constraints!(constraintSet::LinearConstraintSet,permutation::Array{Int,1})::Nothing
    constraintSet.A = constraintSet.A[permutation,:]
    constraintSet.loBs = constraintSet.loBs[permutation]
    constraintSet.upBs = constraintSet.upBs[permutation]
    return
end
