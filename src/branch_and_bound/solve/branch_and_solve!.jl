# @Author: Massimo De Mauri <massimo>
# @Date:   2019-06-04T16:11:53+02:00
# @Email:  massimo.demauri@gmail.com
# @Filename: branch_and_solve!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-05T18:27:42+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}

function branch_and_solve!(node::BBnode,workspace::BBworkspace{T1,T2})::Array{BBnode,1} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # create a list of children
    if node.avgFrac != 0.0
        children = branch!(node,workspace)
    else
        children = [deepcopy(node)]
    end
    # solve all the children
    for k in 1:length(children)
        solve!(children[k],workspace)
        #TODO compute pseudoCost
    end
    # return the solved children
    return children
end


function branch!(node::BBnode,workspace::BBworkspace{T1,T2})::Array{BBnode,1} where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # compute the variables fractionality
    fractionality = @. abs(round(node.primal[workspace.dscIndices]) - node.primal[workspace.dscIndices])
    @. fractionality = fractionality*(fractionality>workspace.settings.integerTolerance)

    # select a branching index
    tmpIndex = workspace.settings.branchingPriorityRule(fractionality,workspace.pseudoCosts)
    branchIndex = workspace.dscIndices[tmpIndex]

    # check if the selected variable belongs to a sos1 group
    if length(workspace.sos1Groups)>0 && workspace.sos1Groups[tmpIndex] != -1

        # collect all the variables belonging to the same sos1 groups
        sos1Group = [i for i in 1:length(workspace.dscIndices)
                          if workspace.sos1Groups[i] == workspace.sos1Groups[tmpIndex] &&
                          node.branchLoBs[i] != node.branchUpBs[i]]

        sos1Branching = true
    else
        sos1Branching = false
    end

    numChildren = 0
    if sos1Branching == true && length(sos1Group) > 1 # SOS1 branching

        # number of children to create
        numChildren = 2

        # order the variables in the sos1 group by their fractionality
        permute!(sos1Group,sortperm(fractionality[sos1Group]))

        # create the new bounds for the children nodes
        newLoBs = [copy(node.branchLoBs),copy(node.branchLoBs)]
        newUpBs = [copy(node.branchUpBs),copy(node.branchUpBs)]
        @. newLoBs[1][sos1Group[1:2:end]] = newUpBs[1][sos1Group[1:2:end]] = 0.
        @. newLoBs[2][sos1Group[2:2:end]] = newUpBs[2][sos1Group[2:2:end]] = 0.


    elseif fractionality[tmpIndex] >= workspace.settings.integerTolerance # standard branching

        # number of children to create
        numChildren = 2

        # create the new bounds for the children nodes
        newLoBs = [copy(node.branchLoBs),copy(node.branchLoBs)]
        newLoBs[1][tmpIndex] = ceil(node.primal[branchIndex]-get_primalTolerance(workspace.subsolverWS))
        newUpBs = [copy(node.branchUpBs),copy(node.branchUpBs)]
        newUpBs[2][tmpIndex] = floor(node.primal[branchIndex]+get_primalTolerance(workspace.subsolverWS))

    end

    # create the list of children
    children = Array{BBnode}(undef,numChildren)
    for k in 1:numChildren

        children[k] = BBnode(newLoBs[k],newUpBs[k],copy(node.primal),
                             copy(node.bndDual),copy(node.cnsDual),
                             NaN,NaN,false)
    end

    return children
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
    (objVal,ssStatus,~) = solve!(workspace.subsolverWS,varLoBs,varUpBs,node.primal,node.bndDual,node.cnsDual)

    if ssStatus == 0
        node.objVal = objVal
        node.reliable = true
    elseif ssStatus == 1
        node.objVal = Inf
        node.reliable = true
    elseif ssStatus == 2
        @warn "Numerical difficulties in the subsolver"
        node.objVal = objVal
        node.reliable = false
    else
        @error "Error in the subsolver"
    end

    # count how many subproblems we have solved
    workspace.status.numExploredNodes = workspace.status.numExploredNodes + 1

    # compute node average fractionality and pseudo_cost
    fractionality = @. node.primal[workspace.dscIndices] - round(node.primal[workspace.dscIndices])
    @. fractionality =  fractionality*(fractionality>workspace.settings.integerTolerance)

    node.avgFrac = 2. * sum(@. abs(fractionality))/length(workspace.dscIndices)
    return
end
