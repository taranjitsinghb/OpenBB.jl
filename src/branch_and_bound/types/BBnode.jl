# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:44+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnode.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-23T19:31:09+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractBBnode end; struct NullBBnode <: AbstractBBnode end

# this is the set of data that distinguishes a node from another
mutable struct BBnode <: AbstractBBnode
    varLoBs::Array{Float64,1}
    varUpBs::Array{Float64,1}
    cnsLoBs::Array{Float64,1}
    cnsUpBs::Array{Float64,1}
    primal::Array{Float64,1}
    bndDual::Array{Float64,1}
    cnsDual::Array{Float64,1}
    avgAbsFrac::Float64
    objective::Float64
    pseudoObjective::Float64
    reliable::Bool
end


# construct a BBnode given its lower bounds, upper bounds and solution hotstart
function BBnode(varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},
                cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                primal::Array{Float64,1},bndDual::Array{Float64,1},cnsDual::Array{Float64,1})::BBnode

    return BBnode(varLoBs,varUpBs,cnsLoBs,cnsUpBs,primal,bndDual,cnsDual,NaN,NaN,NaN,false)
end


# this node is used to arrest the branch and bound process
struct KillerNode <:AbstractBBnode
    count::Int
end

# overload of functions
import Base.copy
function copy(node::BBnode)::BBnode
    return BBnode(node.varLoBs,node.varUpBs,
                  node.cnsLoBs,node.cnsUpBs,
                  node.primal,node.bndDual,node.cnsDual,
                  node.avgAbsFrac,node.objective,node.pseudoObjective,node.reliable)
end


import Base.deepcopy
function deepcopy(node::BBnode)::BBnode
    return BBnode(copy(node.varLoBs),copy(node.varUpBs),
                  copy(node.cnsLoBs),copy(node.cnsUpBs),
                  copy(node.primal),copy(node.bndDual),copy(node.cnsDual),
                  node.avgAbsFrac,node.objective,node.pseudoObjective,node.reliable)
end
