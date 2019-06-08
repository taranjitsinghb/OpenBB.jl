# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-06T14:53:54+02:00
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
function run!(workspace::BBworkspace{T1,T2})::Nothing where T1<:AbstractWorkspace where T2<:AbstractSharedMemory

    # timing
    lastTimeCheckpoint = time()

    # update algorithm status
    workspace.status.description = "running"

    # set the algorithm in "active state"
    idle = false

    # id of the process the workspace is being manipulated from
    processId = myid()

    # print the initial algorithm status
    printCountdown = 0.0

    # keep in memory how many nodes the algorithm has explored since the last time
    # it has sent a node to the neighbouring process
    workBalaceCounter = 0

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

            # print algorithm status
            if workspace.settings.verbose && processId == 1
               if printCountdown <= 0
                   print_status(workspace)
                   printCountdown = workspace.settings.statusInfoPeriod
              end
            end

            # apply rounding heuristics (#TODO  rewrite it!)
            # if node.avgFrac <= workspace.settings.roundingHeuristicsThreshold
            #     push!(workspace.activeQueue,simple_rounding_heuristics(node,workspace))
            # end

            # pick a node to process from the activeQueue
            node = pop!(workspace.activeQueue)

            # check if node is already suboptimal (the upper-bound might have changed)
            if node.objVal > workspace.status.objUpB - workspace.settings.primalTolerance
                if workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes are stored
                    push!(workspace.unactivePool,node)
                end

            elseif !(workspace.sharedMemory isa NullSharedMemory) && # we are multiprocessing
                   workBalaceCounter >= 1 && # time to send the node to another process
                   !isready(workspace.sharedMemory.outputChannel) # send only one node per time

                    # send the new node to the neighbouring process
                    put!(workspace.sharedMemory.outputChannel,node)

                    # reset work balance counter
                    workBalaceCounter = 0

            else

                # update work balance counter
                workBalaceCounter += 1

                # create a list of children nodes
                children = branch_and_solve!(node,workspace)

                # remove the infeasible children
                filter!((child)->child.objVal<Inf,children)

                # handle the rest of the children
                for child in children

                    if workspace.status.objUpB < child.objVal + workspace.settings.primalTolerance # the child is suboptimal
                        if workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes are stored
                            push!(workspace.unactivePool,child)
                        end

                    elseif child.avgFrac == 0.0 # a new solution has been found

                        if child.reliable # the solution is reliable

                            # insert new solution into the solutionPool
                            push!(workspace.solutionPool,child)

                            # update the number of solutions found
                            workspace.status.numSolutions += 1
                            # update the objective upper bound
                            workspace.status.objUpB = child.objVal

                            # update the global objective upper bound and the number of solutions found
                            if !(workspace.sharedMemory isa NullSharedMemory)
                                workspace.sharedMemory.stats[1] +=1
                                if workspace.status.objUpB < workspace.sharedMemory.objectiveBounds[end]
                                    workspace.sharedMemory.objectiveBounds[end] = workspace.status.objUpB
                                end
                            end

                            # recompute optimality gaps
                            if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                                workspace.status.absoluteGap = workspace.status.relativeGap = Inf
                            else
                                workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                                workspace.status.relativeGap = workspace.status.absoluteGap/abs(1e-10 + workspace.status.objUpB)
                            end

                        else # the solution is not reliable

                            # store the unreliable solution in the unactivePool
                            push!(workspace.unactivePool,child)
                        end

                    else  # the child has to be inserted in the active_queue

                        # insert the child in the queue
                        insert_node!(workspace.activeQueue,child,workspace.settings.expansionPriorityRule,
                                     workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)
                     end
                end
            end

            # recompute the objective lower bound
            if workspace.status.objLoB == -Inf || node.objVal == workspace.status.objLoB
                # compute the new lower-bound
                newObjLoB = workspace.status.objUpB
                for i in length(workspace.activeQueue):-1:1
                    # if we have unreliable problems in the activeQueue we cannot update the lower bound
                    if !workspace.activeQueue[i].reliable
                        newObjLoB = workspace.status.objLoB
                        break
                    elseif workspace.activeQueue[i].objVal < newObjLoB
                        newObjLoB = workspace.activeQueue[i].objVal
                    end
                end

                if workspace.status.objLoB > newObjLoB + workspace.settings.primalTolerance
                    @warn "branch and bound: the objective lower bound has decreased from "*string(workspace.status.objLoB)*" to "*string(newObjLoB)*"..."
                end

                # update the objective lower-bound
                workspace.status.objLoB = newObjLoB

                # recompute optimality gaps
                if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                    workspace.status.absoluteGap = workspace.status.relativeGap = Inf
                else
                    workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                    workspace.status.relativeGap = workspace.status.absoluteGap/abs(1e-10 + workspace.status.objUpB)
                end
            end
        end

        # communicate with the other processes
        if !(workspace.sharedMemory isa NullSharedMemory)

            if idle
                # declare the worker ready to stop
                workspace.sharedMemory.arrestable[processId] = true

                # wait for a message to arrive in the inputChannel
                while !isready(workspace.sharedMemory.inputChannel) &&
                      !all(@. workspace.sharedMemory.arrestable)
                    # pause for some time
                    pause(1e-5)
                end
            end

            if isready(workspace.sharedMemory.inputChannel)

                # undeclare the worker ready to stop
                if idle
                    workspace.sharedMemory.arrestable[processId] = false
                end

                # go active
                idle = false

                # take a new node from the input channel
                newNode = take!(workspace.sharedMemory.inputChannel)

                # insert the new node in the queue
                insert_node!(workspace.activeQueue,newNode,workspace.settings.expansionPriorityRule,
                             workspace.status,unreliablePriority=workspace.settings.unreliable_subps_priority)

                # update the objective lower bound
                if newNode.objVal < workspace.status.objLoB
                    workspace.status.objLoB = newNode.objVal
                    # recompute optimality gaps
                    if workspace.status.objUpB == Inf || workspace.status.objLoB == -Inf
                        workspace.status.absoluteGap = workspace.status.relativeGap = Inf
                    else
                        workspace.status.absoluteGap = workspace.status.objUpB - workspace.status.objLoB
                        workspace.status.relativeGap = workspace.status.absoluteGap/abs(1e-10 + workspace.status.objUpB)
                    end
                end
            end


            # receive the number of solutions found and the global upper bound
            workspace.status.objUpB = workspace.sharedMemory.objectiveBounds[end]
            workspace.status.numSolutions = workspace.sharedMemory.stats[1]

            # comunicate the current local lower bound
            workspace.sharedMemory.objectiveBounds[processId] = workspace.status.objLoB
        end
    end

    ############################## termination ##############################

    # check for erroneous termination (to remove at release)
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

        workspace.status.description = 3
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
