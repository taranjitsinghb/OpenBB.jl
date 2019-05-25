# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-25T18:40:35+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# insert a list of equivalent nodes into a queue
function insert_nodes!(queue::Array{BBnode,1},nodes::Array{BBnode,1},
                       priorityRule::Function,status::BBstatus;
                       unreliablePriority::Int=0)::Nothing

    if nodes[1].reliable || unreliablePriority == 0
    # normal queue insertion
        tmpIndex = 0
        for i in length(queue):-1:1
            if priorityRule(nodes[1],queue[i],status)
                tmpIndex = i
                break
            end
        end
        splice!(queue,tmpIndex+1:tmpIndex,nodes)

    elseif unreliablePriority == -1
        # put the new unreliable nodes at the bottom of the activeQueue to deal with them later
        splice!(workspace.activeQueue,1:0,out[2])

    elseif unreliablePriority == 1
        # put the new unreliable nodes at the top of the activeQueue to try to fastly get rid of them
        append!(workspace.activeQueue,out[2])
    else
        @error "wrong priority setting for unreliable nodes"
    end

    return
end





# branch and bound algorithm
function run!(workspace::BBworkspace)::Nothing

    # timing
    lastTimeCheckpoint = time()

    # update algorithm status
    workspace.status.description = "running"

    # set the algorithm in "active state"
    idle = false

    # it is necessary to update the objective lowerbound
    objLoBMayHaveChanged = false

    # main loop
    while !idle

        # update total time
        workspace.status.totalTime += time() - lastTimeCheckpoint; lastTimeCheckpoint = time()

        # stopping conditions
        @sync if length(workspace.activeQueue) == 0 || # no more nodes in the queue
           workspace.status.totalTime >= workspace.settings.timeLimit || # time is up
           workspace.settings.custom_stopping_rule(workspace) || # custom stopping rule triggered
           workspace.status.absoluteGap <= workspace.settings.absoluteGapTolerance || # reached required absolute gap
           workspace.status.relativeGap <= workspace.settings.relativeGapTolerance || # reached required relative gap
           (workspace.settings.numSolutionsLimit > 0 && workspace.status.numSolutions >= workspace.settings.numSolutionsLimit) # the required number of solutions has been found

           # set the algorithm in idle state
           idle = true

       else # continue with branch and bound

           # update algorithm status and print it
           if workspace.settings.verbose
               if workspace.settings.iterationInfoFreq == 1 ||
                  mod(workspace.status.numRelaxationsSolved,workspace.settings.iterationInfoFreq) == 1
                   print_status(workspace)
              end
           end

           # pick a node to process from the activeQueue
           node = pop!(workspace.activeQueue)

           # check if it is necessary to update the lowerbound
           if node.objVal == workspace.status.objLoB
               objLoBMayHaveChanged = true
           else
               objLoBMayHaveChanged = false
           end

           # pick a node to send to the neighbouring process from the activeQueue
           if !(workspace.sharedMemory isa NullSharedMemory) && # multiprocessing?
              !isready(workspace.sharedMemory.outputChannel) && # send only one node per time
              length(workspace.activeQueue) > 0 && # there is a node to send
              workspace.activeQueue[end].objVal < workspace.status.objUpB # do not send suboptimal nodes

                # enforce recomputation of lowerbound
                if workspace.activeQueue[end].objVal == workspace.status.objLoB
                    objLoBMayHaveChanged = true
                end

                # send a new node to the neighbouring process
                nodeToSend = pop!(workspace.activeQueue)
                @async put!(workspace.sharedMemory.outputChannel,nodeToSend)

           end

            # solve the node
            out = solve_and_branch!(node,workspace)

            if out[1] == "solution" && out[2][1].reliable # a reliable solution has been found
                # insert new solution into the solutionPool
                push!(workspace.solutionPool,out[2][1])

                # update the number of solutions found
                workspace.status.numSolutions += 1
                # update the objective upper bound
                workspace.status.objUpB = out[2][1].objVal

                # update the global objective upper bound and the number of solutions found
                if !(workspace.sharedMemory isa NullSharedMemory) && out[2][1].objVal < workspace.sharedMemory.objectiveBounds[end]
                    workspace.sharedMemory.objectiveBounds[end] = out[2][1].objVal
                    workspace.sharedMemory.stats[1] +=1
                end

            elseif out[1] == "solution"  # a not reliable solution has been found
                # store the obtained (not reliable) solution
                push!(workspace.unactivePool,out[2][1])

            elseif out[1] == "suboptimal" && workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes cannot be completely eliminated
                # store the suboptimal node
                push!(workspace.unactivePool,out[2][1])

            elseif out[1] == "children"  # no solution found. two children nodes have been created

                # insert the new nodes into the activeQueue
                insert_nodes!(workspace.activeQueue,out[2],workspace.settings.expansion_priority_rule,
                              workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)

                # apply rounding heuristics
                if node.avgFrac <= workspace.settings.roundingHeuristicsThreshold
                    heuristicnode = simple_rounding_heuristics(node,workspace)
                    push!(workspace.activeQueue,heuristicnode)
                end
            end

            # recompute the lower bound if:
            # 1 - there is no lower bound
            # 2 - the node providing the lower bound has been removed from the activeQueue
            if workspace.status.objLoB == -Inf || objLoBMayHaveChanged

                newObjLoB = workspace.status.objUpB
                for i in length(workspace.activeQueue):-1:1
                    # if we have unreliable problems in the activeQueue we cannot update the lower bound
                    if !workspace.activeQueue[i].reliable
                        newObjLoB = workspace.status.objLoB
                        break
                    elseif workspace.activeQueue[i].objVal < newObjLoB
                        newObjLoB =  workspace.activeQueue[i].objVal
                    end
                end
                if workspace.status.objLoB > newObjLoB + workspace.settings.primalTolerance
                    println("branch and bound: the objective lower bound has decreased from "*string(workspace.status.objLoB)*" to "*string(newObjLoB)*"...")
                end
                workspace.status.objLoB = newObjLoB
            end
        end

        # communicate with the other processes
        if !(workspace.sharedMemory isa NullSharedMemory)

            # check arrest conditions and start arrest procedure
            if idle && # nothing to do locally
               !isready(workspace.sharedMemory.inputChannel) && # no nodes to pick
               workspace.sharedMemory.stats[2] == workspace.settings.numProcesses-1 # all the other workers are waiting.

                 # send a killer-node to the neighbouring process
                 @async put!(workspace.sharedMemory.outputChannel,KillerNode(1))
            end

            # check if a new node is available
            while idle || isready(workspace.sharedMemory.inputChannel)

                if idle
                    # communicate the waiting state
                    workspace.sharedMemory.stats[2] += 1.
                end
                while !isready(workspace.sharedMemory.inputChannel)
                    sleep(0.001)
                end
                # take a new node from the input channel
                newNode = take!(workspace.sharedMemory.inputChannel)

                if idle
                    # exit waiting state
                    workspace.sharedMemory.stats[2] -= 1.
                end

                if newNode isa KillerNode # handle killer-nodes

                    if (idle && newNode.count < workspace.settings.numProcesses-1) ||
                       newNode.count < 2*(workspace.settings.numProcesses-1)
                       while isready(workspace.sharedMemory.inputChannel)
                           sleep(0.001)
                       end
                        # propagate the killerNode
                        @async put!(workspace.sharedMemory.outputChannel,KillerNode(newNode.count+1))
                    end

                    if newNode.count >= workspace.settings.numProcesses-1
                        # go idle (to force arrest) and exit the inner loop
                        idle = true; break
                    end

                else # a normal node: insert it in the queue
                    # insert the new node in the queue
                    insert_nodes!(workspace.activeQueue,[newNode],workspace.settings.expansion_priority_rule,
                                  workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)
                    # update the objective lower bound
                    if newNode.objVal < workspace.status.objLoB
                        workspace.status.objLoB = newNode.objVal
                    end
                    # go active (to prevent arrest) and exit the inner loop
                    idle = false; break
                end
            end

            # update the number of solutions found and the upper bound
            workspace.status.objUpB = workspace.sharedMemory.objectiveBounds[end]
            workspace.status.objLoB = min(workspace.status.objLoB,workspace.status.objUpB)
            workspace.status.numSolutions = workspace.sharedMemory.stats[1]
        end

        # recompute optimality gaps
        if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
            workspace.status.absoluteGap = workspace.status.relativeGap = Inf
        else
            workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
            workspace.status.relativeGap = workspace.status.absoluteGap/abs(1e-10 + workspace.status.objUpB)
        end
    end

    ############################## termination ##############################

    if !(workspace.sharedMemory isa NullSharedMemory)
        # empty the communicationChannels
        while isready(workspace.sharedMemory.inputChannel)
            @assert take!(workspace.sharedMemory.inputChannel) isa KillerNode
        end
    end



    if workspace.status.absoluteGap < workspace.settings.absoluteGapTolerance ||
       workspace.status.relativeGap < workspace.settings.relativeGapTolerance

        workspace.status.description = "optimalSolutionFound"
        if workspace.settings.verbose
            print_status(workspace)
            println(" Exit: Optimal Solution Found")
        end

    elseif length(workspace.activeQueue) == 0 && workspace.status.objUpB == Inf

        workspace.status.description = "infeasible"
        if workspace.settings.verbose
            print_status(workspace)
            println(" Exit: infeasibilty detected")
        end
    # elseif length(workspace.activeQueue) == 0 && length(workspace.solutionPool) > 0
    #
    #     workspace.status.description = "noReliableSolutionFound"
    #     if workspace.settings.verbose
    #         print_status(workspace)
    #         println(" Exit: no reliable solution found")
    #     end

    else
        workspace.status.description = "interrupted"
        if workspace.settings.verbose
            print_status(workspace)
            println(" Exit: interrupted")
        end
    end

    return
end
