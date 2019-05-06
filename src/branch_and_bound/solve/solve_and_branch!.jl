# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T16:19:44+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve_and_branch.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-03T13:49:11+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# this function is called everytime a subproblem is picked from the activeQueue
function solve_and_branch!(subproblem::BBsubproblem, workspace::BBworkspace)::Tuple{String,Array{BBsubproblem,1}}

    # set the subproblem bounds (considering that the general bounds might have changed)
    globalLoBs, globalUpBs = get_variableBounds(workspace)
    varLoBs = copy(globalLoBs)
    varUpBs = copy(globalUpBs)
    for k in keys(subproblem.branchLoBs)
        if varUpBs[k] < subproblem.branchLoBs[k]
            # declare the problem infeasible
            return ("infeasible",Array{BBsubproblem,1}())
        else
            varLoBs[k] = max(varLoBs[k],subproblem.branchLoBs[k])
        end
    end
    for k in keys(subproblem.branchUpBs)
        if varLoBs[k] > subproblem.branchUpBs[k]
            # declare the problem infeasible
            return ("infeasible",Array{BBsubproblem,1}())
        else
            varUpBs[k] = min(varUpBs[k],subproblem.branchUpBs[k])
        end
    end

    for i in 1:length(workspace.dscIndices)
        if varLoBs[workspace.dscIndices[i]] != varUpBs[workspace.dscIndices[i]] && subproblem.pseudoCosts[i] == 0.
            println(subproblem.branchLoBs)
            println(subproblem.branchUpBs)
            println(subproblem.pseudoCosts)
            cacca
        end
    end


    # solve the subproblem
    # subproblem status guide:
    # 0 -> solved
    # 1 -> infeasible
    # 2 -> unreliable
    # 3 -> error
    out = solve!(workspace.subsolverWS;
                    varLoBs=varLoBs,varUpBs=varUpBs,
                    primal=subproblem.primal,
                    bndDual=subproblem.bndDual,cnsDual=subproblem.cnsDual)

    # count how many subproblems we have solved
    workspace.status.numRelaxationsSolved = workspace.status.numRelaxationsSolved + 1

    # check feasibility of the subproblem solution
    if out.status == 1
        return ("infeasible",Array{BBsubproblem,1}())
    end

    # compute subproblem fractionality
    dscv_vars_val = out.primal[workspace.dscIndices]
    fractionality = Array{Float64,1}(undef,length(workspace.dscIndices))
    @. fractionality =  abs(dscv_vars_val - round(dscv_vars_val))
    @. fractionality =  fractionality*(fractionality>workspace.settings.integerTolerance)

    # check the problem for suboptimality, giving some slack to badly solved problems
    if out.objVal + (1 - 101*(out.status==2))*get_primalTolerance(workspace.subsolverWS) > min(workspace.status.objUpB,workspace.settings.objectiveCutoff)

        return ("suboptimal",[BBsubproblem(copy(subproblem.branchUpBs),
                                           copy(subproblem.branchUpBs),
                                           copy(subproblem.pseudoCosts),
                                           copy(out.primal),
                                           copy(out.bndDual),
                                           copy(out.cnsDual),
                                           2*sum(fractionality)/length(workspace.dscIndices),
                                           out.objVal,
                                           out.status== 0)])

    elseif maximum(fractionality) <= workspace.settings.integerTolerance
        # round the integer variables
        @. out.primal[workspace.dscIndices] = round(out.primal[workspace.dscIndices])

        # new solution found!
        return ("solution",[BBsubproblem(copy(subproblem.branchLoBs),
                                         copy(subproblem.branchUpBs),
                                         copy(subproblem.pseudoCosts),
                                         copy(out.primal),
                                         copy(out.bndDual),
                                         copy(out.cnsDual),
                                         2*sum(fractionality)/length(workspace.dscIndices),
                                         out.objVal,
                                         out.status== 0)])
    else

        # select a branching index
        tmpIndex = workspace.settings.branching_priority_rule(fractionality,subproblem.pseudoCosts)
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

            # create the new bounds for the children subproblems
            newLoBs = newUpBs = [Dict{Int,Float64}([workspace.dscIndices[sos1Group[i]] => 0. for i in 1:2:length(sos1Group)]),
                                 Dict{Int,Float64}([workspace.dscIndices[sos1Group[i]] => 0. for i in 2:2:length(sos1Group)])]

            # compute children pseudoCosts (incomplete)
            tmpPseudoCosts = copy(subproblem.pseudoCosts)
            for i in 1:2:length(sos1Group)
                tmpPseudoCosts[sos1Group[i]] = 0.
            end
            childrenPseudoCosts = [tmpPseudoCosts]
            tmpPseudoCosts = copy(subproblem.pseudoCosts)
            for i in 2:2:length(sos1Group)
                tmpPseudoCosts[sos1Group[i]] = 0.
            end
            push!(childrenPseudoCosts,tmpPseudoCosts)

            # compute average fractionality of the children
            childrenAvgFrac = [2*(sum(fractionality)-sum(fractionality[sos1Group[1:2:end]]))/length(workspace.dscIndices),
                               2*(sum(fractionality)-sum(fractionality[sos1Group[2:2:end]]))/length(workspace.dscIndices)]

        elseif fractionality[tmpIndex] >= workspace.settings.integerTolerance # standard branching

            # number of children to create
            numChildren = 2

            # create the new bounds for the children subproblems
            newLoBs = [Dict{Int,Float64}(branchIndex=>ceil(out.primal[branchIndex]-get_primalTolerance(workspace.subsolverWS))),
                       Dict{Int,Float64}()]
            newUpBs = [Dict{Int,Float64}(),
                       Dict{Int,Float64}(branchIndex=>floor(out.primal[branchIndex]+get_primalTolerance(workspace.subsolverWS)))]


            # compute children pseudoCosts (incomplete)
            tmpPseudoCosts = copy(subproblem.pseudoCosts)
            if varUpBs[branchIndex] == newLoBs[1][branchIndex]
                tmpPseudoCosts[tmpIndex] = 0.
            end
            childrenPseudoCosts = [tmpPseudoCosts]
            tmpPseudoCosts = copy(subproblem.pseudoCosts)
            if varLoBs[branchIndex] == newUpBs[2][branchIndex]
                tmpPseudoCosts[tmpIndex] = 0.
            end
            push!(childrenPseudoCosts,tmpPseudoCosts)

            # compute average fractionality of the children
            childrenAvgFrac = repeat([2*(sum(fractionality)-fractionality[tmpIndex])/length(workspace.dscIndices)],2)

        else # case to consider due to possible infinite pseudo costs
        
            # number of children to create (depends on the value of the selected variable)
            numChildren = 1

            # first child: fix the value of the branching variable
            newLoBs = newUpBs = [Dict{Int,Float64}(branchIndex=>round(out.primal[branchIndex]))]
            tmpPseudoCosts = copy(subproblem.pseudoCosts)
            tmpPseudoCosts[tmpIndex] = 0.
            childrenPseudoCosts = [tmpPseudoCosts]


            # second child: round down
            roundedVal = round(out.primal[branchIndex])-1
            if roundedVal >= varLoBs[branchIndex]
                numChildren += 1
                newLoBs = vcat(newLoBs,[Dict{Int,Float64}()])
                newUpBs = vcat(newUpBs,[Dict{Int,Float64}(branchIndex=>roundedVal)])

                tmpPseudoCosts = copy(subproblem.pseudoCosts)
                if roundedVal == varLoBs[branchIndex]
                    tmpPseudoCosts[tmpIndex] = 0.
                end
                push!(childrenPseudoCosts,tmpPseudoCosts)
            end

            # second third child: round up
            roundedVal = round(out.primal[branchIndex])+1
            if roundedVal <= varUpBs[branchIndex]
                numChildren += 1
                newLoBs = vcat(newLoBs,[Dict{Int,Float64}(branchIndex=>roundedVal)])
                newUpBs = vcat(newUpBs,[Dict{Int,Float64}()])

                tmpPseudoCosts = copy(subproblem.pseudoCosts)
                if roundedVal == varUpBs[branchIndex]
                    tmpPseudoCosts[tmpIndex] = 0.
                end
                push!(childrenPseudoCosts,tmpPseudoCosts)
            end

            childrenAvgFrac = repeat([2*(sum(fractionality)-fractionality[tmpIndex])/length(workspace.dscIndices)],numChildren)
        end

        # println(newLoBs," - ",newUpBs)
        # println(subproblem.pseudoCosts)
        # println(out.primal[workspace.dscIndices])

        # create the list of children
        children = Array{BBsubproblem}(undef,numChildren)
        for k in 1:numChildren
            # println(subproblem.branchLoBs,newLoBs[k])
            # println(subproblem.branchUpBs,newUpBs[k])
            # println(childrenPseudoCosts[k])

            # perform bound propagation (TO DO)

            children[k] = BBsubproblem(merge(subproblem.branchLoBs,newLoBs[k]),
                                       merge(subproblem.branchUpBs,newUpBs[k]),
                                       childrenPseudoCosts[k],
                                       copy(out.primal),
                                       copy(out.bndDual),
                                       copy(out.cnsDual),
                                       childrenAvgFrac[k],
                                       out.objVal,
                                       out.status== 0)
        end



        return ("children",children)
    end
end
