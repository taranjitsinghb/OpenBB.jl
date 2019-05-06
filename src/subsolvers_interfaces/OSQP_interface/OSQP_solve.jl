# @Author: Massimo De Mauri <massimodemauri>
# @Date:   2019-02-06T15:17:59+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: OSQP_interface.jl
# @Last modified by:   massimo
# @Last modified time: 2019-03-11T12:25:21+01:00
# @License: apache 2.0
# @Copyright: {{copyright}}


function solve!(workspace::OSQPworkspace;
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
    OSQP.update!(workspace.model;l=vcat(varLoBs,workspace.cnsLoBs),
                                 u=vcat(varUpBs,workspace.cnsUpBs))
    # set hotstart info
    if length(primal) > 0 && length(bndDual) > 0 && length(cnsDual) > 0
        OSQP.warm_start!(workspace.model; x=primal, y=vcat(bndDual,cnsDual))
    end

    # solve problem
    sol = OSQP.solve!(workspace.model)

    # output sol info
    if  sol.info.status_val == 1
        status = 0 # "solved"
    elseif sol.info.status_val == -3
        status = 1 # "infeasible"
    elseif sol.info.status_val in [2,3,4,-6,-2]
        status = 2 # "unreliable"
        sol.x = @. min(max(sol.x,varLoBs),varUpBs)
        @warn "Inaccuracy in subproblem sol, message: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    elseif sol.info.status_val in [-7,-10]
        status = 3 # "error"
        @error "Subsover error, status: "*string(sol.info.status)*" (code: "*string(sol.info.status_val)*")"
    else
        @error "Subsolver unknown status: "*string(sol.info.status)*"("*string(sol.info.status_val)*")"
    end

    #return solution
    nVars = size(workspace.A,2 )
    return SubSolution(sol.x, sol.y[1:nVars],sol.y[nVars+1:end], sol.info.obj_val, status, sol.info.run_time)
end
