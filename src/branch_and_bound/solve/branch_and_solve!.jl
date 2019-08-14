# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-04T16:11:53+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: branch_and_solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-13T20:23:58+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}

function branch_and_solve!(node::BBnode,workspace::BBworkspace{T1,T2})::Array{BBnode,1} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # create a list of children
    if node.avgAbsFrac == 0.0 || isnan(node.objective)
        children, branchIndices_dsc = [deepcopy(node)], [0]
    else
        children, branchIndices_dsc = branch!(node,workspace)
    end


    # solve all the children
    for k in 1:length(children)
        solve!(children[k],workspace)

        # update pseudoCosts
        if branchIndices_dsc[k]>0 && children[k].reliable && children[k].objective < Inf


            # compute objective and primal variation
            deltaObjective = max(children[k].objective-node.objective,workspace.settings.primalTolerance) # the max filters out small numerical errors
            deltaVariable = children[k].primal[workspace.dscIndices[branchIndices_dsc[k]]] - node.primal[workspace.dscIndices[branchIndices_dsc[k]]]

            if deltaVariable < -workspace.settings.primalTolerance

                # update the pseudoCost
                mu = 1/(workspace.pseudoCosts[2][branchIndices_dsc[k],1]+1)
                workspace.pseudoCosts[1][branchIndices_dsc[k],1] = (1-mu)*workspace.pseudoCosts[1][branchIndices_dsc[k],1] - mu*deltaObjective/deltaVariable
                workspace.pseudoCosts[2][branchIndices_dsc[k],1] += 1

            elseif deltaVariable > workspace.settings.primalTolerance

                # update the pseudoCost
                mu = 1/(workspace.pseudoCosts[2][branchIndices_dsc[k],2]+1)
                workspace.pseudoCosts[1][branchIndices_dsc[k],2] = (1-mu)*workspace.pseudoCosts[1][branchIndices_dsc[k],2] + mu*deltaObjective/deltaVariable
                workspace.pseudoCosts[2][branchIndices_dsc[k],2] += 1
            end
        end

    end
    # return the solved children
    return children
end


function branch!(node::BBnode,workspace::BBworkspace{T1,T2})::Tuple{Array{BBnode,1},Array{Int,1}} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # select a branching index
    branchIndex_dsc, priorityScores = branching_priority_rule(workspace.settings.branchingPriorityRule,
                                                       node.primal[workspace.dscIndices],workspace.pseudoCosts,
                                                       workspace.settings.primalTolerance)
    @assert branchIndex_dsc > 0

    # get the index of the branching variable wrt all the variables
    branchIndex = workspace.dscIndices[branchIndex_dsc]

    # check if the selected variable belongs to a sos1 group
    sos1Branching = false
    if length(workspace.sos1Groups)>0 && workspace.sos1Groups[branchIndex_dsc] != -1

        # collect all the variables belonging to the same sos1 groups
        sos1Group = [i for i in 1:length(workspace.dscIndices)
                          if workspace.sos1Groups[i] == workspace.sos1Groups[branchIndex_dsc] &&
                          node.branchLoBs[workspace.dscIndices[i]] != node.branchUpBs[workspace.dscIndices[i]]]

        sos1Branching = true
    end

    # actually branch
    if sos1Branching == true && length(sos1Group) > 1 # SOS1 branching

        # order the variables in the sos1 group by their priority score
        permute!(sos1Group,sortperm((@. abs(node.primal[workspace.dscIndices[sos1Group]])),rev=true))

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),copy(node.bndDual),copy(node.cnsDual))
        @. children[1].branchLoBs[workspace.dscIndices[sos1Group[1:2:end]]] = 0.
        @. children[1].branchUpBs[workspace.dscIndices[sos1Group[1:2:end]]] = 0.

        # second child
        children[2] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),copy(node.bndDual),copy(node.cnsDual))
        @. children[2].branchLoBs[workspace.dscIndices[sos1Group[2:2:end]]] = 0.
        @. children[2].branchUpBs[workspace.dscIndices[sos1Group[2:2:end]]] = 0.

        if length(sos1Group) == 3
            return children, [0,sos1Group[2]]
        elseif length(sos1Group) == 2
            return children, [sos1Group[1],sos1Group[2]]
        else
            return children, [0,0]
        end

    else # standard branching

        # create a list of children nodes
        children = Array{BBnode}(undef,2)

        # first child
        children[1] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),copy(node.bndDual),copy(node.cnsDual))
        children[1].branchLoBs[branchIndex] = ceil(node.primal[branchIndex]-get_primalTolerance(workspace.subsolverWS))

        # second child
        children[2] = BBnode(copy(node.branchLoBs),copy(node.branchUpBs),copy(node.primal),copy(node.bndDual),copy(node.cnsDual))
        children[2].branchUpBs[branchIndex] = floor(node.primal[branchIndex]+get_primalTolerance(workspace.subsolverWS))

        return children, [branchIndex_dsc,branchIndex_dsc]

    end
end


function solve!(node::BBnode,workspace::BBworkspace{T1,T2})::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory


    # set the node bounds (considering that the general bounds might have changed)
    globalLoBs, globalUpBs = get_variableBounds(workspace)
    varLoBs = copy(globalLoBs)
    varUpBs = copy(globalUpBs)
    @. varLoBs = max(varLoBs,node.branchLoBs)
    @. varUpBs = min(varUpBs,node.branchUpBs)

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
        node.pseudoObjective = objective
        for (k,i) in enumerate(workspace.dscIndices)
            node.pseudoObjective +=  min(workspace.pseudoCosts[1][k,1]*(node.primal[i]-floor(node.primal[i]+workspace.settings.primalTolerance)),
                                         workspace.pseudoCosts[1][k,2]*(ceil(node.primal[i]-workspace.settings.primalTolerance)-node.primal[i]))/length(workspace.dscIndices)
        end
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
    workspace.status.numExploredNodes = workspace.status.numExploredNodes + 1

    # compute node average fractionality and pseudo_cost
    absoluteFractionality = @. abs(node.primal[workspace.dscIndices] - round(node.primal[workspace.dscIndices]))
    @. absoluteFractionality =  absoluteFractionality*(absoluteFractionality>workspace.settings.primalTolerance)
    node.avgAbsFrac = 2. * sum(absoluteFractionality)/length(workspace.dscIndices)
    return
end
