# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-07T18:58:51+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: pseudo_costs_initialization.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-11T16:59:21+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}

# wrappers for initialization functions
function initialize_pseudoCosts!(functionTuple::Tuple,pseudoCosts::Array{Float64,2},refNode::BBnode)::Nothing
    @assert functionTuple[1] isa Function
    return functionTuple[1](pseudoCosts,refNode,functionTuple[2:end]...)
end


# actual initialization functions
function initialize_to_constant!(pseudoCosts::Array{Float64,2},refNode::BBnode,value::Float64)::Nothing
    @. pseudoCosts[:,1:2] = value
    return
end

function initialize_with_root_objective!(pseudoCosts::Array{Float64,2},refNode::BBnode,proportionalityConstant::Float64)::Nothing
    @. pseudoCosts[:,1:2]  = max(proportionalityConstant*abs(refNode.objective),1e-4)
    return
end

function initialize_with_strong_branching!(pseudoCosts::Array{Float64,2},refNode::BBnode)::Nothing
    #TODO
    @error "Not Implemented Yet"
    return
end
