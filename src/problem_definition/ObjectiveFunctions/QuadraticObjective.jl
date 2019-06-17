# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:34:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QuadraticObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-23T18:54:29+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# quadratic objective
mutable struct QuadraticObjective{T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}} <: AbstractObjectiveFunction
    Q::T
    L::Array{Float64,1}
end

function QuadraticObjective(;Q::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},L::Array{Float64,1})::QuadraticObjective
    return QuadraticObjective(Q,L)
end

import SparseArrays.sparse
function sparse(objectiveFunction::QuadraticObjective)::QuadraticObjective
    return QuadraticObjective(sparse(objectiveFunction.Q),objectiveFunction.L)
end

function get_sparsity(objectiveFunction::QuadraticObjective{SparseMatrixCSC{Float64,Int}})::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(objectiveFunction)[1:2]
end
