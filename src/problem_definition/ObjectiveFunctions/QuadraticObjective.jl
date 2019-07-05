# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:34:36+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QuadraticObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-02T17:14:20+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}



# quadratic objective
mutable struct QuadraticObjective{T<:Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}}} <: AbstractObjective
    Q::T
    L::Array{Float64,1}
end

function QuadraticObjective(;Q::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},L::Array{Float64,1})::QuadraticObjective
    return QuadraticObjective(Q,L)
end

# ...
function get_numVariables(objective::QuadraticObjective)::Int
    return size(objective.Q,1)
end


# ...
import SparseArrays.sparse
function sparse(objectivection::QuadraticObjective)::QuadraticObjective
    return QuadraticObjective(sparse(objectivection.Q),objectivection.L)
end

# ...
function get_sparsity(objective::QuadraticObjective{SparseMatrixCSC{Float64,Int}})::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(objectivection)[1:2]
end


# ...
function remove_variables!(objective::QuadraticObjective,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objective)))
    objective.Q = objective.Q[toKeep,toKeep]
    objective.L = objective.L[toKeep]
    return
end


# ...
function append_term!(objective::QuadraticObjective,newObjectiveTerm::T)::Nothing where T <: Union{NullObjective,LinearObjective,QuadraticObjective}

    if newObjectiveTerm isa LinearObjective

        numNewVariables = get_numVariables(newObjectiveTerm)
        numOldVariables = get_numVariables(objective)

        objective.Q = vcat(hcat(objective.Q,                            zeros(numOldVariables,numNewVariables)),
                           hcat(zeros(numNewVariables,numOldVariables), zeros(numNewVariables,numNewVariables)))
        objective.L = vcat(objective.L,newObjectiveTerm.L)

    elseif newObjectiveTerm isa QuadraticObjective

        numNewVariables = get_numVariables(newObjectiveTerm)
        numOldVariables = get_numVariables(objective)

        objective.Q = vcat(hcat(objective.Q,                            zeros(numOldVariables,numNewVariables)),
                           hcat(zeros(numNewVariables,numOldVariables), newObjectiveTerm.Q                   ))
        objective.L = vcat(objective.L,newObjectiveTerm.L)

    end

    return
end
