# @Author: Massimo De Mauri <massimo>
# @Date:   2019-05-23T11:26:28+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: BBworkspace.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T13:38:58+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# define abstract and null types
abstract type AbstractWorkspace end; struct NullWorkspace <: AbstractWorkspace end


# this structure collects all the informations needed to run B&B
struct BBworkspace{T1<:AbstractWorkspace,T2<:AbstractSharedMemory} <: AbstractWorkspace
    # problem description
    subsolverWS::T1
    dscIndices::Array{Int64,1}
    sos1Groups::Array{Int64,1}
    pseudoCosts::Array{Float64,1}
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
