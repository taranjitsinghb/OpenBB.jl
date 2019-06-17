# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-04T16:11:53+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: branch_and_solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-17T16:47:44+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function branch_and_solve!(node::BBnode,workspace::BBworkspace{T1,T2})::Array{BBnode,1} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # create a list of children
    if node.avgAbsFrac != 0.0
        children, branchIndices = branch!(node,workspace)
    else
        children, branchIndices = [deepcopy(node)], [0]
    end


    # solve all the children
    for k in 1:length(children)
        solve!(children[k],workspace)

        # update pseudoCosts
        if branchIndices[k]>0 && children[k].reliable && children[k].objective < Inf


            # compute objective and primal variation
            deltaObjective = children[k].objective-node.objective
            deltaVariable = children[k].primal[workspace.dscIndices[branchIndices[k]]] - node.primal[workspace.dscIndices[branchIndices[k]]]

            if deltaVariable < -workspace.settings.primalTolerance

                # update the pseudoCost
                mu = 1/(workspace.pseudoCosts[2][branchIndices[k],1]+1)
                workspace.pseudoCosts[1][branchIndices[k],1] = (1-mu)*workspace.pseudoCosts[1][branchIndices[k],1] - mu*deltaObjective/deltaVariable
                workspace.pseudoCosts[2][branchIndices[k],1] += 1

            elseif deltaVariable > workspace.settings.primalTolerance

                # update the pseudoCost
                mu = 1/(workspace.pseudoCosts[2][branchIndices[k],2]+1)
                workspace.pseudoCosts[1][branchIndices[k],2] = (1-mu)*workspace.pseudoCosts[1][branchIndices[k],2] + mu*deltaObjective/deltaVariable
                workspace.pseudoCosts[2][branchIndices[k],2] += 1
            end
        end

    end
    # return the solved children
    return children
end


function branch!(node::BBnode,workspace::BBworkspace{T1,T2})::Tuple{Array{BBnode,1},Array{Int,1}} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # select a branching index
    tmpIndex, priorityScores = branching_priority_rule(workspace.settings.branchingPriorityRule,
                                                       node.primal[workspace.dscIndices],workspace.pseudoCosts,
                                                       workspace.settings.primalTolerance)
    @assert tmpIndex > 0; branchIndex = workspace.dscIndices[tmpIndex]



    # check if the selected variable belongs to a sos1 group
    sos1Branching = false
    if length(workspace.sos1Groups)>0 && workspace.sos1Groups[tmpIndex] != -1

        # collect all the variables belonging to the same sos1 groups
        sos1Group = [i for i in 1:length(workspace.dscIndices)
                          if workspace.sos1Groups[i] == workspace.sos1Groups[tmpIndex] &&
                          node.branchLoBs[i] != node.branchUpBs[i]]

        sos1Branching = true
    end

    # actually branch
    if sos1Branching == true && length(sos1Group) > 1 # SOS1 branching

        # order the variables in the sos1 group by their priority score
        permute!(sos1Group,sortperm(priorityScores[sos1Group]))

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),
                             copy(node.bndDual),copy(node.cnsDual),NaN,NaN,NaN,false)
        @. children[1].branchLoBs[sos1Group[2:2:end]] = children[1].branchUpBs[sos1Group[2:2:end]] = 0.

        # second child
        children[2] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),
                             copy(node.bndDual),copy(node.cnsDual),NaN,NaN,NaN,false)
        @. children[2].branchLoBs[sos1Group[1:2:end]] = children[2].branchUpBs[sos1Group[1:2:end]] = 0.


        return children, [tmpIndex,tmpIndex] # [sos1Group[2],sos1Group[1]]

    else # standard branching

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),
                             copy(node.bndDual),copy(node.cnsDual),NaN,NaN,NaN,false)
        children[1].branchLoBs[tmpIndex] = ceil(node.primal[branchIndex]-get_primalTolerance(workspace.subsolverWS))

        # second child
        children[2] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),
                             copy(node.bndDual),copy(node.cnsDual),NaN,NaN,NaN,false)
        children[2].branchUpBs[tmpIndex] = floor(node.primal[branchIndex]+get_primalTolerance(workspace.subsolverWS))

        return children, [tmpIndex,tmpIndex]

    end
end


function solve!(node::BBnode,workspace::BBworkspace{T1,T2})::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory


    # set the node bounds (considering that the general bounds might have changed)
    globalLoBs, globalUpBs = get_variableBounds(workspace)
    varLoBs = copy(globalLoBs)
    varUpBs = copy(globalUpBs)
    @. varLoBs[workspace.dscIndices] = max(varLoBs[workspace.dscIndices],node.branchLoBs)
    @. varUpBs[workspace.dscIndices] = min(varUpBs[workspace.dscIndices],node.branchUpBs)

    # solve the node
    # node status guide:
    # 0 -> solved
    # 1 -> infeasible
    # 2 -> unreliable
    # 3 -> error
    if any(@. varLoBs > varUpBs + workspace.settings.primalTolerance)
        ssStatus = 1
    else
        (objective,ssStatus,~) = solve!(workspace.subsolverWS,varLoBs,varUpBs,node.primal,node.bndDual,node.cnsDual)
    end


    if ssStatus == 0
        node.objective = objective
        node.pseudoObjective = objective +
                          maximum(@. (node.primal[workspace.dscIndices]-floor(node.primal[workspace.dscIndices]+workspace.settings.primalTolerance))*workspace.pseudoCosts[1][:,1] +
                                     (ceil(node.primal[workspace.dscIndices]-workspace.settings.primalTolerance)-node.primal[workspace.dscIndices])*workspace.pseudoCosts[1][:,2])
        node.reliable = true
    elseif ssStatus == 1
        node.objective = Inf
        node.pseudoObjective = Inf
        node.reliable = true
    elseif ssStatus == 2
        @warn "Numerical difficulties in the subsolver"
        node.objective = objective
        node.reliable = false
    else
        @error "Error in the subsolver"
    end

    # count how many subproblems we have solved
    workspace.status.numRelaxationsSolved = workspace.status.numRelaxationsSolved + 1

    # compute node average fractionality and pseudo_cost
    fractionality = @. node.primal[workspace.dscIndices] - round(node.primal[workspace.dscIndices])
    @. fractionality =  fractionality*(fractionality>workspace.settings.primalTolerance)

    node.avgAbsFrac = 2. * sum(@. abs(fractionality))/length(workspace.dscIndices)
    return
end
