# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:25:57+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearConstraintSet.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-27T18:14:20+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# constructirs and copy functions (Fundamental. Those are used in Branch and Bound)
# named constructor
function SDPConstraintSet(;A::T,loBs::Array{Float64,2},upBs::Array{Float64,2})::SDPConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return SDPConstraintSet{T}(vec(A),vec(loBs),vec(upBs))
end


import Base.copy
function copy(constraintSet::SDPConstraintSet)::SDPConstraintSet
    return SDPConstraintSet(constraintSet.A)
end

import Base.deepcopy
function deepcopy(constraintSet::SDPConstraintSet)::SDPConstraintSet
    return SDPConstraintSet(deepcopy(constraintSet.A))
end

# type conversion
function SDPConstraintSet{T}(constraintSet::SDPConstraintSet{T})::SDPConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return constraintSet
end

function SDPConstraintSet{T}(constraintSet::SDPConstraintSet)::SDPConstraintSet{T} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return SDPConstraintSet(T(constraintSet.A),constraintSet.loBs,constraintSet.upBs)
end

import SparseArrays.sparse
function sparse(constraintSet::SDPConstraintSet{T})::SDPConstraintSet{SparseMatrixCSC{Float64,Int}} where T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return SDPConstraintSet{SparseMatrixCSC{Float64,Int}}(sparse(constraintSet.A),constraintSet.loBs,constraintSet.upBs)
end


# inspect functions (Fundamental. Those are used in Branch and Bound)
function get_numVariables(constraintSet::SDPConstraintSet)::Int
    return size(constraintSet.A,1)*size(constraintSet.A,2)
end

function get_size(constraintSet::SDPConstraintSet)::Int
    return size(constraintSet.A,1)
end

function get_bounds(constraintSet::SDPConstraintSet)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (constraintSet.loBs,constraintSet.upBs)
end


function get_sparsity(constraintSet::SDPConstraintSet)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(constraintSet.A))[1:2]
end

function get_sparsity(constraintSet::SDPConstraintSet,index::Int)::Array{Int,1}
    return findnz(sparse(constraintSet.A[index,:]))[1]
end

function get_firstNZs(constraintSet::SDPConstraintSet)::Array{Int,1}
    out = Array{Int,1}(undef,size(constraintSet.A,1))
    for k in 1:length(out)
        try
            out[k] = findfirst(!iszero,constraintSet.A[k,:])
        catch
            out[k] = -1
        end
    end
    return out
end

function get_firstNZ(constraintSet::SDPConstraintSet,index::Int)::Int
    return findfirst(!iszero,constraintSet.A[index,:])
end

function get_lastNZs(constraintSet::SDPConstraintSet)::Array{Int,1}
    out = Array{Int,1}(undef,size(constraintSet.A,1))
    for k in 1:length(out)
        try
            out[k] = findlast(!iszero,constraintSet.A[k,:])
        catch
            out[k] = -1
        end
    end
    return out
end

function get_lastNZ(constraintSet::SDPConstraintSet,index::Int)::Int
    return findlast(!iszero,constraintSet.A[index,:])
end

function get_linearConstraints(constraintSet::SDPConstraintSet)::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}
    return constraintSet.A
end

# update functions (Not fundamental. Those are used only in updating the problem)
function update_bounds!(constraintSet::SDPConstraintSet;loBs::Array{Float64,1}=Float64[],upBs::Array{Float64,1}=Float64[])::Nothing
    if length(loBs) > 0
        @assert length(loBs) == length(constraintSet.loBs) == length(constraintSet.upBs)
        constraintSet.loBs = loBs
    end
    if length(upBs) > 0
        @assert length(upBs) == length(constraintSet.loBs) == length(constraintSet.upBs)
        constraintSet.upBs = upBs
    end
    return
end

function remove_variables!(constraintSet::SDPConstraintSet,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(constraintSet)))
    constraintSet.A = constraintSet.A[:,toKeep]
    return
end

function insert_variables!(constraintSet::SDPConstraintSet,numVariables::Int,insertionPoint::Int)::Nothing
    @assert numVariables>=0
    @assert 0<=insertionPoint<=get_numVariables(constraintSet)+1
    constraintSet.A = hcat(constraintSet.A[:,1:insertionPoint-1],zeros(size(constraintSet.A,1),numVariables),constraintSet.A[:,insertionPoint:end])
    return
end

function append_variables!(constraintSet::SDPConstraintSet,numVariables::Int)::Nothing
    insert_variables!(constraintSet,numVariables,get_numVariables(constraintSet)+1)
    return
end

function remove_constraints!(constraintSet::SDPConstraintSet,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_size(constraintSet)))
    constraintSet.A = constraintSet.A[toKeep,:]
    constraintSet.loBs = constraintSet.loBs[toKeep]
    constraintSet.upBs = constraintSet.upBs[toKeep]
    return
end

function insert_constraints!(constraintSet1::SDPConstraintSet,constraintSet2::SDPConstraintSet,insertionPoint::Int)::Nothing
        constraintSet1.A = vcat(constraintSet1.A[1:insertionPoint-1,:],constraintSet2.A,constraintSet1.A[insertionPoint:end,:])
        splice!(constraintSet1.loBs,insertionPoint:insertionPoint-1,copy(constraintSet2.loBs))
        splice!(constraintSet1.upBs,insertionPoint:insertionPoint-1,copy(constraintSet2.upBs))
    return
end

function append_constraints!(constraintSet1::SDPConstraintSet,constraintSet2::SDPConstraintSet)::Nothing
    insert_constraints!(constraintSet1,constraintSet2,get_size(constraintSet1)+1)
    return
end

function permute_constraints!(constraintSet::SDPConstraintSet,permutation::Array{Int,1})::Nothing
    constraintSet.A = constraintSet.A[permutation,:]
    constraintSet.loBs = constraintSet.loBs[permutation]
    constraintSet.upBs = constraintSet.upBs[permutation]
    return
end
