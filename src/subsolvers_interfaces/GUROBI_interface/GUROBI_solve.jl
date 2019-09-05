# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-04T12:30:24+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: GUROBIsolve.jl
# @Last modified by:   massimo
# @Last modified time: 2019-09-03T14:54:30+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

function solve!(node::BBnode,workspace::GUROBIworkspace)::Tuple{Int8,Float64}

    # update the gurobi model
    nVars = get_numVariables(workspace)
    nCnss = get_numConstraints(workspace)

    # create a Gurobi environment
    model = Gurobi.gurobi_model(workspace.environment,H = workspace.Q,
                                                      f = workspace.L,
                                                      A = vcat(-workspace.A,workspace.A),
                                                      b = vcat(-node.cnsLoBs,node.cnsUpBs),
                                                      lb = node.varLoBs,
                                                      ub = node.varUpBs)

    Gurobi.update_model!(model)

    # solve problem
    runtime = @elapsed Gurobi.optimize(model)
    status = Gurobi.get_status_code(model)


    # output sol info
    if  status == 2
        status = 0 # "solved"
        node.primal = Gurobi.get_solution(model)
        node.bndDual = zeros(nVars)
        node.cnsDual = zeros(nCnss)
        node.objGap = workspace.settings.OptimalityTol
        node.objVal = (transpose(node.primal)*workspace.Q*node.primal)/2. + transpose(workspace.L)*node.primal
    elseif status in [3,4]
        status = 1 # "infeasible"
        node.primal = NaNs(nVars)
        node.bndDual = NaNs(nVars)
        node.cnsDual = NaNs(nCnss)
        node.objGap = workspace.settings.OptimalityTol
        node.objVal = Inf
    elseif status in [7,8,10,11,13]
        status = 2 # "unreliable"
        node.primal = Gurobi.get_solution(model)
        node.primal = @. min(max(node.primal,node.varLoBs),node.varUpBs)
        node.bndDual = zeros(nVars)
        node.cnsDual = zeros(nCnss)
        newObjVal = (transpose(node.primal)*workspace.Q*node.primal)/2. + transpose(workspace.L)*node.primal
        if newObjVal >= node.ObjVal - node.objGap
            node.objGap = newObjVal - node.objVal + node.objGap
            node.objVal = newObjVal
        else
            node.objGap = Inf # gurobi doesn't give enough information to estimate the gap
            node.objVal = newObjVal
            @warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(status)*")"
        end
    elseif status in [1,5,12]
        status = 3 # "error"
        @error "Subsover error, status: "*string(Gurobi.get_status(model))*" (code: "*string(status)*")"
    else
        @error "Subsolver unknown status: "*string(Gurobi.get_status(model))*" (code:"*string(status)*")"
    end

    return (status, runtime)
end
