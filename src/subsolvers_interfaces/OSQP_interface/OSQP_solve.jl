# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T15:17:59+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-02T12:28:13+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function solve!(subsolverWS::OSQPworkspace,
                varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                primal::Array{Float64,1},bndDual::Array{Float64,1},cnsDual::Array{Float64,1})::Tuple{Float64,Int8,Float64}

    # collect info on the problem
    numVars = get_numVariables(subsolverWS)

    # update the bounds in the osqp model
    OSQP.update!(subsolverWS.model;l=vcat(varLoBs,cnsLoBs),u=vcat(varUpBs,cnsUpBs))

    # set hotstart info
    if length(primal) > 0 && length(bndDual) > 0 && length(cnsDual) > 0
        OSQP.warm_start!(subsolverWS.model; x=copy(primal), y=vcat(bndDual,cnsDual))
    end

    # solve problem
    sol = OSQP.solve!(subsolverWS.model)

    # output sol info
    if  sol.info.status_val == 1
        status = 0 # "solved"
        obj_val = sol.info.obj_val
        @. primal = sol.x
        @. bndDual = sol.y[1:numVars]
        @. cnsDual = sol.y[numVars+1:end]
    elseif sol.info.status_val == -3
        status = 1 # "infeasible"
        obj_val = Inf
        @. primal = @. min(max(sol.x,varLoBs),varUpBs)
        @. bndDual = sol.y[1:numVars]
        @. cnsDual = sol.y[numVars+1:end]
    elseif sol.info.status_val in [2,4,3,-6,-2]
        status = 2 # "unreliable
        @. primal = min(max(sol.x,varLoBs),varUpBs)
        @. bndDual = sol.y[1:numVars]
        @. cnsDual = sol.y[numVars+1:end]
        obj_val = 1/2 * transpose(primal) * subsolverWS.Q *primal + transpose(subsolverWS.L) * primal
        @warn "Inaccuracy in node sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    return (obj_val, status, sol.info.run_time)
end
