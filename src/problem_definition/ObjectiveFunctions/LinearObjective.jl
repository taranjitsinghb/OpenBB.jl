# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T18:33:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: LinearObjective.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-05T15:11:30+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# linear objective
mutable struct LinearObjective <: AbstractObjective
    L::Array{Float64,1}
end

function LinearObjective(;L::Array{Float64,1})::LinearObjective
    return LinearObjective(L)
end


# ...
function get_numVariables(objectiveFun::LinearObjective)::Int
    return size(objectiveFun.L,1)
end


# ...
function remove_variables!(objective::LinearObjective,indices::Array{Int,1})::Nothing
    toKeep = filter(x->!(x in indices), collect(1:get_numVariables(objective)))
    objective.L = objective.L[toKeep]
    return
end


# ...
function append_term!(objective::LinearObjective,newObjectiveTerm::T)::Nothing where T <: Union{NullObjective,LinearObjective}
    if newObjectiveTerm isa LinearObjective
        objective.L = vcat(objective.L,newObjectiveTerm.L)
    end
    return
end
