# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-26T15:03:06+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: ObjectiveFunctions.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-17T22:59:50+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# abstract and null types
abstract type AbstractObjectiveFun end; struct NullObjectiveFun <: AbstractObjectiveFun end

# linear objective
mutable struct LinearObj <: AbstractObjectiveFun
    L::Array{Float64,1}
end

function LinearObj(;L::Array{Float64,1})::LinearObj
    return LinearObj(L)
end

# quadratic objective
mutable struct QuadraticObj{T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}} <: AbstractObjectiveFun
    Q::T
    L::Array{Float64,1}
end

function QuadraticObj(;Q::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},L::Array{Float64,1})::QuadraticObj
    return QuadraticObj(Q,L)
end

import Base.size
function get_numVariables(objectiveFunction::Union{LinearObj,QuadraticObj})::Int
    return length(objectiveFunction.L)
end

import SparseArrays.sparse
function sparse(objectiveFunction::QuadraticObj)::QuadraticObj
    return QuadraticObj(sparse(objectiveFunction.Q),objectiveFunction.L)
end

function get_sparsity(objectiveFunction::QuadraticObj{SparseMatrixCSC{Float64,Int}})::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(objectiveFunction)[1:2]
end
