# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T15:17:59+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QPALM_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-24T14:43:24+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function solve!(subsolverWS::QPALMworkspace,
                varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                primal::Array{Float64,1},bndDual::Array{Float64,1},cnsDual::Array{Float64,1})::Tuple{Float64,Int8,Float64}

    # update the qpalm model
    QPALM.update!(subsolverWS.model;bmin=vcat(varLoBs,cnsLoBs),bmax=vcat(varUpBs,cnsUpBs))

    # set hotstart info
    if length(primal) > 0 && length(bndDual) > 0 && length(cnsDual) > 0
        QPALM.warm_start!(subsolverWS.model; x_warm_start=primal, y_warm_start=vcat(bndDual,cnsDual))
    end

    # solve problem
    sol = QPALM.solve!(subsolverWS.model)

    # compute the objective value
    if sol.info.status_val in [1,2,3,4,-6,-2]
        obj_val = 1/2 * transpose(sol.x) * subsolverWS.Q * sol.x + transpose(subsolverWS.L) * sol.x
    else
        obj_val = NaN
    end

    # output sol info
    if  sol.info.status_val == 1
        status = 0 # "solved"
    elseif sol.info.status_val == -3
        status = 1 # "infeasible"
    elseif sol.info.status_val in [2,3,4,-6,-2]
        status = 2 # "unreliable"
        sol.x = @. min(max(sol.x,varLoBs),varUpBs)
        @warn "Inaccuracy in node sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    #return solution
    nVars = get_numVariables(subsolverWS)
    @. primal = sol.x
    @. bndDual = sol.y[1:nVars]
    @. cnsDual = sol.y[nVars+1:end]

    return (obj_val, status, sol.info.run_time)
end
