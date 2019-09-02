# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T15:17:59+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T14:56:40+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function solve!(node::BBnode,workspace::OSQPworkspace)::Tuple{Int8,Float64}

    # collect info on the problem
    numVars = get_numVariables(workspace)

    # update the bounds in the osqp model
    OSQP.update!(workspace.model;l=vcat(node.varLoBs,node.cnsLoBs),u=vcat(node.varUpBs,node.cnsUpBs))

    # set hotstart info
    if length(node.primal) > 0 && length(node.bndDual) > 0 && length(node.cnsDual) > 0
        OSQP.warm_start!(workspace.model; x=copy(node.primal), y=vcat(node.bndDual,node.cnsDual))
    end

    # solve problem
    sol = OSQP.solve!(workspace.model)

    # output sol info
    if  sol.info.status_val == 1
        status = 0 # "solved"
        node.objective = sol.info.obj_val
        @. node.primal = sol.x
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:end]
    elseif sol.info.status_val == -3
        status = 1 # "infeasible"
        node.objective = Inf
        @. node.primal = @. min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:end]
    elseif sol.info.status_val in [2,4,3,-6,-2]
        status = 2 # "unreliable
        @. node.primal = min(max(sol.x,node.varLoBs),node.varUpBs)
        @. node.bndDual = sol.y[1:numVars]
        @. node.cnsDual = sol.y[numVars+1:end]
        node.objective = 1/2 * transpose(node.primal) * workspace.Q *node.primal + transpose(workspace.L) * node.primal
        @warn "Inaccuracy in node sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    return (status, sol.info.run_time)
end
