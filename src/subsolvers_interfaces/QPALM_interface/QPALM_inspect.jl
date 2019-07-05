# @Author: Massimo De Mauri <massimo>
# @Date:   2019-03-16T15:53:30+01:00
# @Email:  massimo.demauri@gmail.com
# @Filename: QPALM_inspect.jl
# @Last modified by:   massimo
# @Last modified time: 2019-07-01T13:04:40+02:00
# @License: LGPL-3.0
# @Copyright: {{copyright}}

# ...
function get_solver_name(workspace::QPALMworkspace)::String
    return "QPALM"
end

# this function returns the sparsity pattern of the constraint set
function get_constraints_sparsity(workspace::QPALMworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(workspace.A))[1:2]
end

# this function returns the sparsity pattern a constraint in the constraint set
function get_constraint_sparsity(workspace::QPALMworkspace,index::Int)::Array{Int,1}
    return findnz(sparse(workspace.A[index,:]))[1]
end

# this function returns the sparsity pattern of the objevtive function
function get_objective_sparsity(workspace::QPALMworkspace)::Tuple{Array{Int,1},Array{Int,1}}
    return findnz(sparse(workspace.Q))[1:2]
end

#
function get_numVariables(workspace::QPALMworkspace)::Int
    return size(workspace.A,2)
end

#
function get_numConstraints(workspace::QPALMworkspace)::Int
    return size(workspace.A,1)
end


#
function get_variableBounds(workspace::QPALMworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (workspace.varLoBs,workspace.varUpBs)
end

#
function get_constraintBounds(workspace::QPALMworkspace)::Tuple{Array{Float64,1},Array{Float64,1}}
    return (workspace.cnsLoBs,workspace.cnsUpBs)
end

# ...
function get_settings(workspace::QPALMworkspace)::QPALMsettings
    return workspace.settings
end

# ...
function get_constraints(workspace::QPALMworkspace)::LinearConstraintSet
    return LinearConstraintSet(A=workspace.A,loBs=workspace.cnsLoBs,upBs=workspace.cnsUpBs)
end

# ...
function get_objective(workspace::QPALMworkspace)::QuadraticObjective
    return QuadraticObjective(Q=workspace.Q,L=workspace.L)
end
