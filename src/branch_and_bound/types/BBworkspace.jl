# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:28+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBworkspace.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T14:49:57+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractWorkspace end; struct NullWorkspace <: AbstractWorkspace end


# this structure collects all the informations needed to run B&B
mutable struct BBworkspace{T1<:AbstractWorkspace,T2<:AbstractSharedMemory} <: AbstractWorkspace
    # problem description
    subsolverWS::T1
    dscIndices::Array{Int64,1}
    sos1Groups::Array{Int64,1}
    pseudoCosts::Tuple{Array{Float64,2},Array{Int,2}}
    # branch and bound status
    activeQueue::Array{BBnode,1}
    solutionPool::Array{BBnode,1}
    unactivePool::Array{BBnode,1}
    status::BBstatus
    # multiprocessing communication
    sharedMemory::T2
    # user settings
    settings::BBsettings
end



# construct the root node for a problem of the given dimensions
function BBroot(workspace::BBworkspace)::BBnode
    numVars = get_numVariables(workspace)
    numCnss = get_numConstraints(workspace)
    varBounds = get_variableBounds(workspace)
    cnsBounds = get_constraintBounds(workspace)

    return BBnode(copy(varBounds[1]),copy(varBounds[2]),
                  copy(cnsBounds[1]),copy(cnsBounds[2]),
                  zeros(numVars),zeros(numVars),zeros(numCnss))
end

# construct the root node for a problem of the given dimensions (with primal initialization)
function BBroot(workspace::BBworkspace,primal::Array{Float64,1})::BBnode
    numVars = get_numVariables(workspace)
    numCnss = get_numConstraints(workspace)
    varBounds = get_variableBounds(workspace)
    cnsBounds = get_constraintBounds(workspace)

    return BBnode(copy(varBounds[1]),copy(varBounds[2]),
                  copy(cnsBounds[1]),copy(cnsBounds[2]),
                  copy(primal),zeros(numVars),zeros(numCnss))
end
