# @Author: Massimo De Mauri <massimo>
# @Date:   2019-02-19T20:30:02+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: types.jl
# @Last modified by:   massimo
# @Last modified time: 2019-05-06T17:31:03+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}


# structure used to hold the settings for QPALM
mutable struct QPALMsettings <: AbstractSettings
    max_iter::Int                       # Maximum number of interations 0 < max_iter
    eps_abs::Float64                    # Absolute tolerance 	0 <= eps_abs
    eps_rel::Float64                    # Relative tolerance 	0 <= eps_rel
    eps_abs_in::Float64                 # Intermediate absolute convergence tolerance 	0 <= eps_abs_in
    eps_rel_in::Float64                 # Intermediate relative convergence tolerance 	0 <= eps_rel_in
    rho::Float64                        # Tolerance scaling factor step 	0 < rho < 1
    eps_prim_inf::Float64               # Primal infeasibility tolerance    0 <= eps_prim_inf
    eps_dual_inf::Float64               # Dual infeasibility tolerance      0 <= eps_dual_inf
    theta::Float64                      # Penalty update criterion          theta <= 1
    delta::Float64                      # Penalty update factor             1 < delta
    tau_init::Float64                   # Initial stepsize in backtracking  0 < tau_init
    proximal::Bool                      # Use proximal method of multipliers  True/False
    gamma_init::Float64                 # Injtial proximal penalty          0 < gamma_init
    gamma_upd::Float64                  # Proximal penalty update factor    1 <= gamma_upd
    gamma_max::Float64                  # Proximal penalty parameter cap    gamma_init < gamma_max
    scaling::Int                        # Number of scaling iterations 	0 (disabled) or 0 < scaling (integer)
    nonconvex::Bool                     # Indicate if QP is convex          True/False
    verbose::Bool                       # Print output 	True/False
    warm_start::Bool                    # Perform warm starting 	True/False
end


function QPALMsettings(; max_iter::Int=10000,
                         eps_abs::Float64=1e-4,
                         eps_rel::Float64=1e-4,
                         eps_abs_in::Float64=1.0,
                         eps_rel_in::Float64=1.0,
                         rho::Float64=0.1,
                         eps_prim_inf::Float64=1e-4,
                         eps_dual_inf::Float64=1e-4,
                         theta::Float64=0.25,
                         delta::Float64=10.0,
                         tau_init::Float64=1.0,
                         proximal::Bool=true,
                         gamma_init::Float64=1e6,
                         gamma_upd::Float64=10.0,
                         gamma_max::Float64=1e8,
                         scaling::Int=10,
                         nonconvex::Bool=false,
                         verbose::Bool=false,
                         warm_start::Bool=true)::QPALMsettings



    return QPALMsettings(max_iter,eps_abs,eps_rel,eps_abs_in,eps_rel_in,rho,eps_prim_inf,
                        eps_dual_inf,theta,delta,tau_init,proximal,gamma_init,
                        gamma_upd,gamma_max,scaling,nonconvex,verbose,warm_start)
end





# structure used for storing data for QPALM solver
mutable struct QPALMworkspace <: AbstractWorkspace
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
    model::QPALM.Model
    settings::QPALMsettings
end


# getters and setters
function get_primalTolerance(workspace::QPALMworkspace)::Float64
    return workspace.settings.eps_prim_inf
end
