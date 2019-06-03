# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:44+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve_and_branch.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-03T16:36:14+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# this function is called everytime a node is picked from the activeQueue
function solve_and_branch!(node::BBnode, workspace::BBworkspace)::Tuple{String,Array{BBnode,1}}

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
    out = solve!(workspace.subsolverWS;
                    varLoBs=varLoBs,varUpBs=varUpBs,
                    primal=node.primal,bndDual=node.bndDual,cnsDual=node.cnsDual)

    # count how many nodes we have solved
    workspace.status.numExploredNodes = workspace.status.numExploredNodes + 1
    # check feasibility of the node solution
    if out.status == 1
        return ("infeasible",Array{BBnode,1}())
    end

    # compute node fractionality
    dscv_vars_val = out.primal[workspace.dscIndices]
    fractionality = Array{Float64,1}(undef,length(workspace.dscIndices))
    @. fractionality =  abs(dscv_vars_val - round(dscv_vars_val))
    @. fractionality =  fractionality*(fractionality>workspace.settings.integerTolerance)

    # check the problem for suboptimality, giving some slack to badly solved problems
    if out.objVal + (1 - 101*(out.status==2))*get_primalTolerance(workspace.subsolverWS) > min(workspace.status.objUpB,workspace.settings.objectiveCutoff)

        return ("suboptimal",[BBnode(copy(node.branchLoBs),copy(node.branchUpBs),
                                     copy(out.primal),copy(out.bndDual),copy(out.cnsDual),
                                     2*sum(fractionality)/length(workspace.dscIndices),
                                     out.objVal,out.status==0)])

    elseif maximum(fractionality) <= workspace.settings.integerTolerance
        # round the integer variables
        @. out.primal[workspace.dscIndices] = round(out.primal[workspace.dscIndices])

        # new solution found!
        return ("solution",[BBnode(copy(node.branchLoBs),copy(node.branchUpBs),
                                   copy(out.primal),copy(out.bndDual),copy(out.cnsDual),
                                   2*sum(fractionality)/length(workspace.dscIndices),
                                   out.objVal,out.status==0)])
    else

        # select a branching index
        tmpIndex = workspace.settings.branching_priority_rule(fractionality,workspace.pseudoCosts)
        branchIndex = workspace.dscIndices[tmpIndex]


        # check if the selected variable belongs to a sos1 group
        if length(workspace.sos1Groups)>0 && workspace.sos1Groups[tmpIndex] != -1

            # collect all the variables belonging to the same sos1 groups
            sos1Group = [i for i in 1:length(workspace.dscIndices)
                            if workspace.sos1Groups[i] == workspace.sos1Groups[tmpIndex] &&
                               varLoBs[workspace.dscIndices[i]] != varUpBs[workspace.dscIndices[i]]]

            sos1Branching = true
        else
            sos1Branching = false
        end


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

            # compute average fractionality of the children
            childrenAvgFrac = [2*(sum(fractionality)-sum(fractionality[sos1Group[1:2:end]]))/length(workspace.dscIndices),
                               2*(sum(fractionality)-sum(fractionality[sos1Group[2:2:end]]))/length(workspace.dscIndices)]

        elseif fractionality[tmpIndex] >= workspace.settings.integerTolerance # standard branching

            # number of children to create
            numChildren = 2

            # create the new bounds for the children nodes
            newLoBs = [copy(node.branchLoBs),copy(node.branchLoBs)]
            newLoBs[1][tmpIndex] = ceil(out.primal[branchIndex]-get_primalTolerance(workspace.subsolverWS))
            newUpBs = [copy(node.branchUpBs),copy(node.branchUpBs)]
            newUpBs[2][tmpIndex] = floor(out.primal[branchIndex]+get_primalTolerance(workspace.subsolverWS))

            # compute average fractionality of the children
            childrenAvgFrac = repeat([2*(sum(fractionality)-fractionality[tmpIndex])/length(workspace.dscIndices)],2)

        end

        # create the list of children
        children = Array{BBnode}(undef,numChildren)
        for k in 1:numChildren

            # perform bound propagation (TO DO)

            children[k] = BBnode(newLoBs[k],newUpBs[k],copy(out.primal),
                                 copy(out.bndDual),copy(out.cnsDual),
                                 childrenAvgFrac[k],out.objVal,out.status==0)
        end



        return ("children",children)
    end
end
