# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-06T18:33:16+01:00
# @Email:  massimo.demauri@gmail.com
# @Project: OpenBB
# @Filename: run!.jl
# @Last modified by:   massimo
# @Last modified time: 2019-06-17T14:24:43+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# insert a list of equivalent nodes into a queue
function insert_node!(queue::Array{BBnode,1},node::BBnode,
                       priorityRule::Tuple,status::BBstatus;
                       unreliablePriority::Int=0)::Nothing

    if node.reliable || unreliablePriority == 0
    # normal queue insertion
        insertionPoint = 1
        for i in length(queue):-1:1
            if expansion_priority_rule(priorityRule,node,queue[i],status)
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
    timeToShareNodes = false

    # main loop
    while !idle

        # get elapsed time for the last iteration
        iterationTime = time() - lastTimeCheckpoint
        lastTimeCheckpoint = time()
        # update total time
        workspace.status.totalTime += iterationTime
        printCountdown -= iterationTime

        # print algorithm status
        if workspace.settings.verbose && processId == 1
           if printCountdown <= 0
               print_status(workspace)
               printCountdown = workspace.settings.statusInfoPeriod
          end
        end


        # stopping conditions
        if length(workspace.activeQueue) == 0 || # no more nodes in the queue
           workspace.status.totalTime >= workspace.settings.timeLimit || # time is up
           workspace.settings.customStoppingRule(workspace) || # custom stopping rule triggered
           workspace.status.absoluteGap <= workspace.settings.absoluteGapTolerance || # reached required absolute gap
           workspace.status.relativeGap <= workspace.settings.relativeGapTolerance || # reached required relative gap
           (workspace.settings.numSolutionsLimit > 0 && workspace.status.numSolutions >= workspace.settings.numSolutionsLimit) # the required number of solutions has been found

           # set the algorithm in idle state
           idle = true

       else # continue with branch and bound

            # apply rounding heuristics (#TODO  rewrite it!)
            # if node.avgAbsFrac <= workspace.settings.roundingHeuristicsThreshold
            #     push!(workspace.activeQueue,simple_rounding_heuristics(node,workspace))
            # end

            # pick a node to process from the activeQueue
            node = pop!(workspace.activeQueue)

            # check if node is already suboptimal (the upper-bound might have changed)
            if node.objective > workspace.status.objUpB - workspace.settings.primalTolerance
                if workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes are stored
                    push!(workspace.unactivePool,node)
                end

            elseif !(workspace.sharedMemory isa NullSharedMemory) && # we are multiprocessing
                   timeToShareNodes # it is time to send a child to the next process
                   # !isready(workspace.sharedMemory.outputChannel) # the channel is free

                # send the new node to the neighbouring process
                put!(workspace.sharedMemory.outputChannel,node)

                # the next node has to be explored locally
                timeToShareNodes = false

            else # explore the node locally

                # the next node will be sent to the neighbouring process
                timeToShareNodes = true

                # create a list of children nodes
                children = branch_and_solve!(node,workspace)

                # remove the infeasible children
                filter!((child)->child.objective<Inf,children)

                # handle the rest of the children
                for child in children

                    if workspace.status.objUpB < child.objective + workspace.settings.primalTolerance # the child is suboptimal
                        if workspace.settings.dynamicMode # in dynamic mode the suboptimal nodes are stored
                            push!(workspace.unactivePool,child)
                        end

                    elseif child.avgAbsFrac == 0.0 # a new solution has been found

                        if child.reliable # the solution is reliable

                            # insert new solution into the solutionPool
                            push!(workspace.solutionPool,child)

                            # update the number of solutions found
                            workspace.status.numSolutions += 1
                            # update the objective upper bound
                            workspace.status.objUpB = child.objective

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
                    else # insert the child in the local queue
                        # insert the child in the queue
                        insert_node!(workspace.activeQueue,child,workspace.settings.expansionPriorityRule,
                                     workspace.status,unreliablePriority=workspace.settings.unreliableSubproblemsPriority)
                    end
                end
            end

            # recompute the objective lower bound
            if workspace.status.objLoB == -Inf || node.objective == workspace.status.objLoB
                # compute the new lower-bound
                newObjLoB = workspace.status.objUpB
                for i in length(workspace.activeQueue):-1:1
                    # if we have unreliable problems in the activeQueue we cannot update the lower bound
                    if !workspace.activeQueue[i].reliable
                        newObjLoB = workspace.status.objLoB
                        break
                    elseif workspace.activeQueue[i].objective < newObjLoB
                        newObjLoB = workspace.activeQueue[i].objective
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

                # keep track of the time spent in waiting
                waitingStartTime = time()

                # wait for a message to arrive in the inputChannel
                while !isready(workspace.sharedMemory.inputChannel) &&
                      !all(@. workspace.sharedMemory.arrestable)
                    # pause for some time
                    pause(1e-3)
                end

                # keep track of the time spent in waiting
                workspace.status.waitingTime += time()-waitingStartTime
            end

            if isready(workspace.sharedMemory.inputChannel)
                if idle
                    # undeclare the worker ready to stop
                    workspace.sharedMemory.arrestable[processId] = false
                    # go active
                    idle = false
                end

                # take a new node from the input channel
                newNode = take!(workspace.sharedMemory.inputChannel)

                # insert the new node in the queue
                insert_node!(workspace.activeQueue,newNode,workspace.settings.expansionPriorityRule,
                             workspace.status,unreliablePriority=workspace.settings.unreliableSubproblemsPriority)

                # update the objective lower bound
                if newNode.objective < workspace.status.objLoB
                    workspace.status.objLoB = newNode.objective
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
            workspace.status.objLoB = min(workspace.status.objLoB,workspace.status.objUpB)
            workspace.status.numSolutions = workspace.sharedMemory.stats[1]

            # comunicate the current local lower bound
            workspace.sharedMemory.objectiveBounds[processId] = workspace.status.objLoB
        end
    end

    ############################## termination ##############################

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
            println(" Exit: Infeasibilty Detected")
        end

    else
        workspace.status.description = "interrupted"
        if workspace.settings.verbose && processId == 1
            print_status(workspace)
            println(" Exit: Interrupted")
        end
    end

    return
end
