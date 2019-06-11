# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:44+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBnode.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-11T13:32:18+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractBBnode end; struct NullBBnode <: AbstractBBnode end

# this is the set of data that distinguishes a node from another
mutable struct BBnode <: AbstractBBnode
    branchLoBs::Array{Float64,1}
    branchUpBs::Array{Float64,1}
    primal::Array{Float64,1}
    bndDual::Array{Float64,1}
    cnsDual::Array{Float64,1}
    avgAbsFrac::Float64
    objective::Float64
    pseudoObjective::Float64
    reliable::Bool
end

# this node is used to arrest the branch and bound process
struct KillerNode <:AbstractBBnode
    count::Int
end

# overload of functions
import Base.copy
function copy(node::BBnode)::BBnode
    return BBnode(node.branchLoBs,node.branchUpBs,
                  node.primal,node.bndDual,node.cnsDual,
                  node.avgAbsFrac,node.objective,node.pseudoObjective,node.reliable)
end


import Base.deepcopy
function deepcopy(node::BBnode)::BBnode
    return BBnode(copy(node.branchLoBs),copy(node.branchUpBs),
                  copy(node.primal),copy(node.bndDual),copy(node.cnsDual),
                  node.avgAbsFrac,node.objective,node.pseudoObjective,node.reliable)
end
