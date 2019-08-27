# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-04T12:30:24+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: GUROBIsolve.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-27T15:36:19+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


function solve!(subsolverWS::GUROBIworkspace,
                varLoBs::Array{Float64,1},varUpBs::Array{Float64,1},cnsLoBs::Array{Float64,1},cnsUpBs::Array{Float64,1},
                primal::Array{Float64,1},bndDual::Array{Float64,1},cnsDual::Array{Float64,1})::Tuple{Float64,Int8,Float64}

    # update the gurobi model
    nVars = get_numVariables(subsolverWS)
    nCnss = get_numConstraints(subsolverWS)

    model = Gurobi.Model(subsolverWS.environment,"OpenBBsubsolver")
    if nVars > 0
        Gurobi.add_cvars!(model,subsolverWS.L,varLoBs,varUpBs)
        Gurobi.update_model!(model)
        Gurobi.add_qpterms!(model,sparse(subsolverWS.Q))
        Gurobi.add_rangeconstrs!(model,sparse(subsolverWS.A),cnsLoBs,cnsUpBs)
        Gurobi.update_model!(model)
    end

    # solve problem
    info_runtime = @elapsed Gurobi.optimize(model)
    info_status = Gurobi.get_status_code(model)


    # output sol info
    if  info_status == 2
        status = 0 # "solved"
        info_primal = Gurobi.get_solution(model)[1:nVars]
        info_bndDual = zeros(nVars)
        info_cnsDual = zeros(nCnss)
        info_obj = (transpose(info_primal)*subsolverWS.Q*info_primal)/2. + transpose(subsolverWS.L)*info_primal

    elseif info_status in [3,4]
        status = 1 # "infeasible"
        info_primal = NaNs(nVars)
        info_bndDual = NaNs(nVars)
        info_cnsDual = NaNs(nCnss)
        info_obj = Inf

    elseif info_status in [7,8,10,11,13]
        status = 2 # "unreliable"
        info_primal = Gurobi.get_solution(model)[1:nVars]
        info_primal = @. min(max(info_primal,varLoBs),varUpBs)
        info_bndDual = zeros(nVars)
        info_cnsDual = zeros(nCnss)
        info_obj = (transpose(info_primal)*subsolverWS.Q*info_primal)/2. + transpose(subsolverWS.L)*info_primal
        @warn "Inaccuracy in node sol, status: "*string(sol.info.status)*" (code: "*string(info_status)*")"

    elseif info_status in [1,5,12]
        status = 3 # "error"
        @error "Subsover error, status: "*string(Gurobi.get_status(model))*" (code: "*string(info_status)*")"
    else
        @error "Subsolver unknown status: "*string(Gurobi.get_status(model))*" (code:"*string(info_status)*")"
    end

    @. primal = info_primal
    @. bndDual = info_bndDual
    @. cnsDual = info_cnsDual
    return (info_obj, status, info_runtime)
end
