# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-31T13:46:31+02:00
# @License: apache 2.0
# @Copyright: {{copyright}}


# insert a list of equivalent nodes into a queue
function insert_node!(queue::Array{BBnode,1},node::BBnode,
                       priorityRule::Function,status::BBstatus;
                       unreliablePriority::Int=0)::Nothing

    if node.reliable || unreliablePriority == 0
    # normal queue insertion
        insertionPoint = 1
        for i in length(queue):-1:1
            if priorityRule(node,queue[i],status)
                insertionPoint = i+1
                break
            end
        end
        insert!(queue,insertionPoint,node)

    elseif unreliablePriority == -1
        # put the new unreliable nodes at the bottom of the activeQueue to deal with them later
        insert!(workspace.activeQueue,0,node)

    elseif unreliablePriority == 1
        # put the new unreliable nodes at the top of the activeQueue to try to fastly get rid of them
        append!(workspace.activeQueue,node)
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

    # id of the process the workspace is being manipulated from
    processId = myid()

    # print the initial algorithm status
    printCountdown = 0.0

    # main loop
    while !idle

        # update total time
        elapsedTime = time() - lastTimeCheckpoint
        workspace.status.totalTime += elapsedTime
        printCountdown -= elapsedTime
        lastTimeCheckpoint = time()

        # stopping conditions
        if length(workspace.activeQueue) == 0 || # no more nodes in the queue
           workspace.status.totalTime >= workspace.settings.timeLimit || # time is up
           workspace.settings.custom_stopping_rule(workspace) || # custom stopping rule triggered
           workspace.status.absoluteGap <= workspace.settings.absoluteGapTolerance || # reached required absolute gap
           workspace.status.relativeGap <= workspace.settings.relativeGapTolerance || # reached required relative gap
           (workspace.settings.numSolutionsLimit > 0 && workspace.status.numSolutions >= workspace.settings.numSolutionsLimit) # the required number of solutions has been found

           # set the algorithm in idle state
           idle = true

       else # continue with branch and bound

           # update print algorithm status
           if workspace.settings.verbose && processId == 1
               if printCountdown <= 0
                   print_status(workspace)
                   printCountdown = workspace.settings.statusInfoPeriod
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

                #  send one of the children to the neighbouring node
                if !(workspace.sharedMemory isa NullSharedMemory) && # multiprocessing?
                   !isready(workspace.sharedMemory.outputChannel) && # send only one node per time
                   length(out[2]) > 1 # there actually is a node to send

                    for k in 1:length(out[2])

                        if k == 2
                           # send a new node to the neighbouring process
                           put!(workspace.sharedMemory.outputChannel,out[2][2])
                       else
                           # insert the remaining children into the queue
                           insert_node!(workspace.activeQueue,out[2][k],workspace.settings.expansion_priority_rule,
                                        workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)
                        end

                    end
               else
                   # insert the new nodes into the activeQueue
                   for k in 1:length(out[2])
                       insert_node!(workspace.activeQueue,out[2][k],workspace.settings.expansion_priority_rule,
                                    workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)
                   end
               end

                # apply rounding heuristics
                if node.avgFrac <= workspace.settings.roundingHeuristicsThreshold
                    push!(workspace.activeQueue,simple_rounding_heuristics(node,workspace))
                end
            end

            # recompute the objective lower bound
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
                 put!(workspace.sharedMemory.outputChannel,KillerNode(1))

            end

            # check if a new node is available
            while idle || isready(workspace.sharedMemory.inputChannel)

                if idle
                    # communicate the waiting state
                    workspace.sharedMemory.stats[2] += 1.
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

                        # propagate the killerNode
                        put!(workspace.sharedMemory.outputChannel,KillerNode(newNode.count+1))
                    end

                    if newNode.count >= workspace.settings.numProcesses-1
                        # go idle (to force arrest) and exit the inner loop
                        idle = true; break
                    end

                else # a normal node: insert it in the queue
                    # insert the new node in the queue
                    insert_node!(workspace.activeQueue,newNode,workspace.settings.expansion_priority_rule,
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
            workspace.status.numSolutions = workspace.sharedMemory.stats[1]

            # comunicate the current lowerbound
            workspace.sharedMemory.objectiveBounds[processId] = workspace.status.objLoB
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
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: Optimal Solution Found")
        end

    elseif length(workspace.activeQueue) == 0 && workspace.status.objUpB == Inf

        workspace.status.description = "infeasible"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: infeasibilty detected")
        end
    # elseif length(workspace.activeQueue) == 0 && length(workspace.solutionPool) > 0
    #
    #     workspace.status.description = "noReliableSolutionFound"
    #     if workspace.settings.verbose && processId == 1
    #         print_status(workspace)
    #         println(" Exit: no reliable solution found")
    #     end

    else
        workspace.status.description = "interrupted"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: interrupted")
        end
    end

    return
end
