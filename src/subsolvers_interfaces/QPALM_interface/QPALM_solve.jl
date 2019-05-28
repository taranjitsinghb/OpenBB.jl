# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T15:17:59+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QPALM_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-11T12:25:21+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function solve!(workspace::QPALMworkspace;
                    varLoBs::Array{Float64,1}=Float64[],
                    varUpBs::Array{Float64,1}=Float64[],
                    primal::Array{Float64,1}=Float64[],
                    bndDual::Array{Float64,1}=Float64[],
                    cnsDual::Array{Float64,1}=Float64[])::SubSolution

    # check inputs
    if length(varLoBs) == 0
        varLoBs = workspace.varLoBs
    end
    if length(varUpBs) == 0
        varUpBs = workspace.varUpBs
    end

    # update the osqp model
    QPALM.update!(workspace.model;bmin=vcat(varLoBs,workspace.cnsLoBs),
                                 bmax=vcat(varUpBs,workspace.cnsUpBs))
    # set hotstart info
    if length(primal) > 0 && length(bndDual) > 0 && length(cnsDual) > 0
        QPALM.warm_start!(workspace.model; x_warm_start=primal, y_warm_start=vcat(bndDual,cnsDual))
    end

    # solve problem
    sol = QPALM.solve!(workspace.model)

    if sol.info.status_val in [1,2,3,4,-6,-2]
        obj_val = 1/2 * transpose(sol.x) * workspace.Q * sol.x + transpose(workspace.L) * sol.x
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
        obj_val = 0 # TODO
        @warn "Inaccuracy in node sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    #return solution
    nVars = size(workspace.A,2 )
    return SubSolution(sol.x, sol.y[1:nVars],sol.y[nVars+1:end], obj_val, status, sol.info.run_time)
end
