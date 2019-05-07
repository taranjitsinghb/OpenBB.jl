# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: solve_singlecore!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-07T15:18:10+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# branch and bound algorithm in single core mode
function run!(workspace::BBworkspace)::Nothing


    # timing
    last_time = time()

    # update algorithm status
    workspace.status.description = "running"

    if false # debug only
        for p in workspace.activeQueue
            println(p.branchLoBs,p.branchUpBs)
        end
    end

    # main loop
    while length(workspace.activeQueue) > 0 &&
          workspace.status.totalTime + time() - last_time <= workspace.settings.timeLimit &&
          !workspace.settings.custom_stopping_rule(workspace) &&
          workspace.status.absoluteGap >= workspace.settings.absoluteGapTolerance &&
          workspace.status.relativeGap >= workspace.settings.relativeGapTolerance &&
          (workspace.settings.numSolutionsLimit == 0 || workspace.status.numSolutions < workspace.settings.numSolutionsLimit)

        # update algorithm status and print it
        workspace.status.totalTime += time() - last_time; last_time = time()
        workspace.status.totalIterations += 1
        if workspace.settings.verbose
            if workspace.settings.iterationInfoFreq == 1 || mod(workspace.status.totalIterations,workspace.settings.iterationInfoFreq) == 1
               print_status(workspace)
           end
        end

        # take a problem from the activeQueue
        subproblem = pop!(workspace.activeQueue)

        if false # debug only
            println("subproblem")
            println("- lBnd: ",subproblem.branchLoBs)
            println("- uBnd: ",subproblem.branchUpBs)
            println("- primal: ",subproblem.primal)
            println("- objV: ",subproblem.objVal," - avgF: ",subproblem.avgFrac)
        end

        # solve the subproblem
        out = solve_and_branch!(subproblem,workspace)


        if false # debug only
            if out[1] == "infeasible"
                println("\n  - infeasible")
            elseif out[1] == "solution"
                println("\n  - solution")
            else
                println("\n  - ",out[1])
                for prb in out[2]
                    # println(" - prim: ",prb.primal[workspace.dscIndices])
                    println("   lBnd:",prb.branchLoBs)
                    println("   uBnd:",prb.branchUpBs)
                    println()
                end
            end
        end


        if out[1] == "solution" && out[2][1].reliable # a reliable solution has been found

            # update the number of solutions found
            workspace.status.numSolutions = workspace.status.numSolutions + 1

            # update the objective upper bound
            if out[2][1].objVal < workspace.settings.objectiveCutoff
                workspace.status.objUpB = out[2][1].objVal
            end

            # insert new solution into the solutionPool
            tmpIndex = 0
            for i in length(workspace.solutionPool):-1:1
                if workspace.settings.expansion_priority_rule(out[2][1],workspace.solutionPool[i],workspace.status)
                    tmpIndex = i
                    break
                end
            end
            splice!(workspace.solutionPool,tmpIndex+1:tmpIndex,out[2])



        elseif out[1] == "solution"  # a not reliable solution has been found
            # store the obtained (not reliable) solution
            push!(workspace.unactivePool,out[2][1])


        elseif out[1] == "suboptimal" && workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes cannot be completely eliminated
            # store the suboptimal subproblem
            push!(workspace.unactivePool,out[2][1])


        elseif out[1] == "children" && out[2][1].reliable # no solution found. two children subproblems have been created

            # insert new problems into the activeQueue
            tmpIndex = 0
            for i in length(workspace.activeQueue):-1:1
                if workspace.settings.expansion_priority_rule(out[2][1],workspace.activeQueue[i],workspace.status)
                    tmpIndex = i
                    break
                end
            end
            splice!(workspace.activeQueue,tmpIndex+1:tmpIndex,out[2])


            # apply rounding heuristics
            if subproblem.avgFrac <= workspace.settings.roundingHeuristicsThreshold
                heuristicSubproblem = simple_rounding_heuristics(subproblem,workspace)
                push!(workspace.activeQueue,heuristicSubproblem)
            end

        elseif out[1] == "children" # no solution found. two not reliable children subproblems have been created

            # deal with the unreliability
            if workspace.settings.unreliable_subps_priority == -1
                # put the new unreliable subproblems at the bottom of the activeQueue to deal with them later
                splice!(workspace.activeQueue,1:0,out[2])

            elseif workspace.settings.unreliable_subps_priority == 0
                # insert the new unreliable subproblems in the activeQueue normally
                tmpIndex = 0
                for i in length(workspace.activeQueue):-1:1
                    if workspace.settings.expansion_priority_rule(out[2][1],workspace.activeQueue[i],workspace.status)
                        tmpIndex = i
                        break
                    end
                end
                splice!(workspace.activeQueue,tmpIndex+1:tmpIndex,out[2])

            elseif workspace.settings.unreliable_subps_priority == 1
                # put the new unreliable subproblems at the top of the activeQueue to try to fastly get rid of them
                append!(workspace.activeQueue,out[2])

            else
                @error "OpenBB, wrong value for \"unreliable_subps_priority\" setting (allowed: {-1,0,1})"
            end
        end


        # recompute the lower bound if:
        # 1 - there is no lower bound
        # 2 - the subproblem providing the lower bound has been removed from the activeQueue
        # 3 - the queue is empty
        if length(workspace.activeQueue) == 0
        # last iteration: only for coherence of results
            workspace.status.objLoB = workspace.status.objUpB

        elseif workspace.status.objLoB == -Inf || subproblem.objVal == workspace.status.objLoB
        # there is no lower bound or the subproblem providing the last lower bound was removed from the queue

            newObjLoB = workspace.status.objUpB
            for i in length(workspace.activeQueue):-1:1
                # if we have not reliable problems in the activeQueue we cannot update the lower bound
                if !workspace.activeQueue[i].reliable
                    newObjLoB = workspace.status.objLoB
                    break
                elseif workspace.activeQueue[i].objVal < newObjLoB
                    newObjLoB =  workspace.activeQueue[i].objVal
                end
            end
            if workspace.status.objLoB > newObjLoB + workspace.settings.primalTolerance
                @error "branch and bound: the objective lower bound has decreased from "*string(workspace.status.objLoB)*" to "*string(newObjLoB)*"... something is wrong"
            end
            workspace.status.objLoB = newObjLoB
        end

        # recompute optimality gaps
        if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
            workspace.status.absoluteGap = workspace.status.relativeGap = Inf
        else
            workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
            workspace.status.relativeGap = workspace.status.absoluteGap/abs(1e-10 + workspace.status.objUpB)
        end

    end

    if length(workspace.activeQueue) == 0 || workspace.status.absoluteGap < workspace.settings.primalTolerance || workspace.status.relativeGap == 0.0

        if length(workspace.solutionPool) > 0
            solution_found = false
            for k in length(workspace.solutionPool):-1:1
                if workspace.solutionPool[k].reliable
                    solution_found = true
                    break
                end
            end

            if solution_found
                workspace.status.description = "best_subproblem_found"
            else
                workspace.status.description = "no_reliable_solution_found"
            end
        else
            workspace.status.description = "no_solution_found"
        end
    else
        workspace.status.description = "interrupted"
    end

    if workspace.settings.verbose print_status(workspace) end

    return
end
