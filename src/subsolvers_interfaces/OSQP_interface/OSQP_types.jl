# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:30:02+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: types.jl
# @Last modified by:   massimo
# @Last modified time: 2019-08-29T17:31:52+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# structure used to hold the settings for OSQP
mutable struct OSQPsettings <: AbstractSettings
    rho::Float64                        # ADMM rho step 	0 < rho
    sigma::Float64                      # ADMM sigma step 	0 < sigma
    max_iter::Int                       # Maximum number of iterations 	0 < max_iter (integer)
    eps_abs::Float64                    # Absolute tolerance 	0 <= eps_abs
    eps_rel::Float64                    # Relative tolerance 	0 <= eps_rel
    eps_prim_inf::Float64               # Primal infeasibility tolerance 	0 <= eps_prim_inf
    eps_dual_inf::Float64               # Dual infeasibility tolerance 	0 <= eps_dual_inf
    alpha::Float64                      # ADMM overrelaxation parameter 	0 < alpha < 2
    delta::Float64                      # Polishing regularization parameter 	0 < delta
    polish::Bool                        # Perform polishing 	True/False
    polish_refine_iter::Int             # Refinement iterations in polish 	0 < polish_refine_iter
    verbose::Bool                       # Print output 	True/False
    scaled_termination::Bool            # Scaled termination conditions 	True/False
    check_termination::Int              # Check termination interval 	0 (disabled) or 0 < check_termination
    warm_start::Bool                    # Perform warm starting 	True/False
    scaling::Int                        # Number of scaling iterations 	0 (disabled) or 0 < scaling (integer)
    adaptive_rho::Bool                  # Adaptive rho 	True/False
    adaptive_rho_interval::Int          # Adaptive rho interval 	0 (automatic) or 0 < adaptive_rho_interval
    adaptive_rho_tolerance::Float64     # Tolerance for adapting rho 	1 <= adaptive_rho_tolerance
    adaptive_rho_fraction::Float64      # Adaptive rho interval as fraction of setup time (auto mode) 	0 < adaptive_rho_fraction
    timeLimit::Float64	                # Run time limit in seconds 	0 (disabled) or 0 <= timeLimit

end


function OSQPsettings(; rho::Float64=1e-1,
                        sigma::Float64=1e-6,
                        max_iter::Int=10000,
                        eps_abs::Float64=1e-6,
                        eps_rel::Float64=1e-6,
                        eps_prim_inf::Float64=1e-4,
                        eps_dual_inf::Float64=1e-4,
                        alpha::Float64=1.6,
                        delta::Float64=1e-06,
                        polish::Bool=true,
                        polish_refine_iter::Int=15,
                        verbose::Bool=false,
                        scaled_termination::Bool=false,
                        check_termination::Int=15,
                        warm_start::Bool=true,
                        scaling::Int=15,
                        adaptive_rho::Bool=true,
                        adaptive_rho_interval::Int=0,
                        adaptive_rho_tolerance::Float64=5.,
                        adaptive_rho_fraction::Float64=0.4,
                        timeLimit::Float64=0.)::OSQPsettings



    return OSQPsettings(rho,sigma,max_iter,eps_abs,eps_rel,eps_prim_inf,eps_dual_inf,alpha,delta,
                        polish,polish_refine_iter,verbose,scaled_termination,check_termination,
                        warm_start,scaling,adaptive_rho,adaptive_rho_interval,adaptive_rho_tolerance,
                        adaptive_rho_fraction,timeLimit)
end





# structure used for storing data for OSQP solver
mutable struct OSQPworkspace <: AbstractWorkspace
    # objective
    Q::SparseMatrixCSC{Float64}
    L::Array{Float64,1}

    # constraints
    A::SparseMatrixCSC{Float64}
    cnsLoBs::Array{Float64,1}
    cnsUpBs::Array{Float64,1}

    # variables
    varLoBs::Array{Float64,1}
    varUpBs::Array{Float64,1}

    # workspace
    model::OSQP.Model
    settings::OSQPsettings
end


# getters and setters
function get_primalTolerance(workspace::OSQPworkspace)::Float64
    return workspace.settings.eps_prim_inf
end
